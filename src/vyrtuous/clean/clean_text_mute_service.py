"""!/bin/python3
clean_text_mute_service.py The purpose of this program is to clean expired text-mutes.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from datetime import datetime, timedelta, timezone

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.text_mute import TextMute

MODEL = TextMute


async def clean_expired_text_mutes() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    expired_text_mutes = await database_factory.select(expired=True, singular=False)
    if expired_text_mutes:
        for expired_text_mute in expired_text_mutes:
            channel_snowflake = int(expired_text_mute.channel_snowflake)
            guild_snowflake = int(expired_text_mute.guild_snowflake)
            member_snowflake = int(expired_text_mute.member_snowflake)
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                bot.logger.debug(
                    f"Unable to locate guild {guild_snowflake} while cleaning up expired text-mute."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                bot.logger.debug(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}) while cleaning up expired text-mute."
                )
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if not simplified_member:
                    bot.logger.debug(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) while cleaning up expired text-mute."
                    )
                    continue
                count += 1
                await database_factory.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                    member_snowflake=member_snowflake,
                )
            else:
                try:
                    await channel.set_permissions(
                        target=member,
                        send_messages=None,
                        add_reactions=None,
                        reason="Cleaning up expired text-mute",
                    )
                except discord.Forbidden as e:
                    bot.logger.error(str(e).capitalize())
                except discord.HTTPException as e:
                    bot.logger.error(f"HTTP error removing expired text mute: {e}")
    return count


async def clean_text_mute_overwrites() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    now = datetime.now(timezone.utc)
    text_mutes = await database_factory.select(singular=False)
    for text_mute in text_mutes:
        channel_snowflake = int(text_mute.channel_snowflake)
        guild_snowflake = int(text_mute.guild_snowflake)
        member_snowflake = int(text_mute.member_snowflake)
        where_kwargs = {
            "channel_snowflake": channel_snowflake,
            "guild_snowflake": guild_snowflake,
            "member_snowflake": member_snowflake,
        }
        set_kwargs = {"reset": True}
        if not text_mute.reset and text_mute.last_muted < now - timedelta(weeks=1):
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                bot.logger.debug(
                    f"Unable to locate guild {guild_snowflake} for removing overwrite."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                bot.logger.debug(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                )
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                bot.logger.debug(
                    f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                )
                continue
            count += 1
            try:
                await channel.set_permissions(
                    target=member,
                    overwrite=None,
                    reason="Resetting text-mute overwrite",
                )
            except discord.Forbidden as e:
                bot.logger.error(str(e).capitalize())
            await database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )
    return count
