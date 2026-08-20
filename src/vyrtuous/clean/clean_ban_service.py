"""!/bin/python3
clean_ban_service.py The purpose of this program is clean expired bans.

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
from vyrtuous.db.ban import Ban
from vyrtuous.db.database_factory import DatabaseFactory

MODEL = Ban


async def clean_expired_bans() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    expired_bans = await database_factory.select(expired=True, singular=False)
    if expired_bans:
        for expired_ban in expired_bans:
            channel_snowflake = int(expired_ban.channel_snowflake)
            guild_snowflake = int(expired_ban.guild_snowflake)
            member_snowflake = int(expired_ban.member_snowflake)
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                bot.logger.debug(
                    f"Unable to locate guild {guild_snowflake} while cleaning up expired ban."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                bot.logger.debug(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake} while cleaning up expired ban."
                )
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if not simplified_member:
                    bot.logger.debug(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) while cleaning up expired ban."
                    )
                    continue
                count += 1
            else:
                try:
                    await channel.set_permissions(
                        member, view_channel=None, reason="Cleaning up expired ban."
                    )
                except discord.Forbidden as e:
                    bot.logger.error(str(e).capitalize())
                except discord.HTTPException as e:
                    bot.logger.error(f"HTTP error removing expired ban: {e}")
            await database_factory.delete(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                member_snowflake=member_snowflake,
            )

    return count


async def clean_ban_overwrites() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    bans = await database_factory.select(singular=False)
    for ban in bans:
        channel_snowflake = int(ban.channel_snowflake)
        guild_snowflake = int(ban.guild_snowflake)
        member_snowflake = int(ban.member_snowflake)
        where_kwargs = {
            "channel_snowflake": channel_snowflake,
            "guild_snowflake": guild_snowflake,
            "member_snowflake": member_snowflake,
        }
        set_kwargs = {"reset": True}
        if (
            not ban.reset
            and ban.last_kicked < datetime.now(timezone.utc) - timedelta(weeks=1)
            and not ban.blacklisted
        ):
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
                    target=member, overwrite=None, reason="Resetting ban overwrite."
                )
            except discord.Forbidden as e:
                bot.logger.error(str(e).capitalize())
            await database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )
    return count
