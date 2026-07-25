"""!/bin/python3
clean_voice_mute_service.py The purpose of this program is to clean expired voice mutes.

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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute

MODEL = VoiceMute


async def clean_expired_voice_mutes() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    expired_voice_mutes = await database_factory.select(expired=True, singular=False)
    if expired_voice_mutes:
        for expired_voice_mute in expired_voice_mutes:
            channel_snowflake = int(expired_voice_mute.channel_snowflake)
            guild_snowflake = int(expired_voice_mute.guild_snowflake)
            member_snowflake = int(expired_voice_mute.member_snowflake)
            target = expired_voice_mute.target
            guild = bot.get_guild(guild_snowflake)
            kwargs = {
                "channel_snowflake": channel_snowflake,
                "guild_snowflake": guild_snowflake,
                "member_snowflake": member_snowflake,
                "target": target,
            }
            if guild is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate guild {guild_snowflake}, cleaning up expired voice-mute."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                )
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.name}), cleaning up expired voice-mute."
                )
                continue
            await database_factory.delete(**kwargs)
            count += 1
            if (
                member.voice
                and member.voice.channel
                and member.voice.channel.id == channel_snowflake
            ):
                try:
                    await member.edit(mute=False)
                    bot.logger.info(
                        f"Undone voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake})."
                    )
                except discord.Forbidden as e:
                    bot.logger.warning(
                        f"Unable to undo voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
                    )
            else:
                bot.logger.info(
                    f"Member {member.display_name} ({member.id}) is not in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), skipping undo voice-mute."
                )
    return count
