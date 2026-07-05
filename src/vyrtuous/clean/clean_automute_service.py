"""!/bin/python3stage"
automute_channel_service.py The purpose of this program is to extend Service to service the stage class.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute

MODEL = AutoMute


async def clean_expired_automutes() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    automute_database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    expired_automutes = await automute_database_factory.select(
        expired=True, singular=False
    )
    voice_mute_database_factory: DatabaseFactory = DatabaseFactory(VoiceMute)
    if expired_automutes:
        for expired_automute in expired_automutes:
            channel_snowflake = int(expired_automute.channel_snowflake)
            guild_snowflake = int(expired_automute.guild_snowflake)
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                await automute_database_factory.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                )
                bot.logger.info(
                    f"Unable to locate guild {guild_snowflake}, cleaning up expired automute."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                await automute_database_factory.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                )
                bot.logger.info(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                )
                continue
            if not isinstance(channel, discord.VoiceChannel):
                continue
            automutes = await voice_mute_database_factory.select(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                target="auto",
                singular=False,
            )
            for automute in automutes:
                member_snowflake = automute.member_snowflake
                member = guild.get_member(member_snowflake)
                if member is None:
                    continue
                else:
                    await voice_mute_database_factory.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                        member_snowflake=member_snowflake,
                        target="auto",
                    )
                    count += 1
                    if member in channel.members:
                        try:
                            await member.edit(mute=False, reason="Undoing automute")
                        except discord.Forbidden:
                            continue
        bot.logger.info("Cleaned up expired automutes.")
    return count
