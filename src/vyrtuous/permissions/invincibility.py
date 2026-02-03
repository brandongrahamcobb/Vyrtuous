"""invincibility.py A utility module for granting invincibility to moderation events and revoking all moderation events for a user.

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

from typing import Dict, Tuple

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.infractions.ban.ban import Ban
from vyrtuous.db.infractions.flag import Flag
from vyrtuous.db.infractions.tmute.text_mute import TextMute
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute
from vyrtuous.utils.logger import logger


class Invincibility:

    state: bool = False
    invincible_members: Dict[Tuple[int, int], bool] = {}

    @classmethod
    async def unrestrict(cls, guild_snowflake, member_snowflake):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "member_snowflake": member_snowflake,
        }
        bans = await Ban.select(**kwargs)
        text_mutes = await TextMute.select(**kwargs)
        voice_mutes = await VoiceMute.select(**kwargs)
        if bans:
            for ban in bans:
                channel = guild.get_channel(ban.channel_snowflake)
                if channel:
                    try:
                        await channel.set_permissions(member, overwrite=None)
                    except discord.Forbidden:
                        logger.warning(
                            f"Unable to unban {member.name} ({member.id}) in {channel.name} ({channel.id})."
                        )
        if text_mutes:
            for text_mute in text_mutes:
                channel = guild.get_channel(text_mute.channel_snowflake)
                if channel:
                    try:
                        await channel.set_permissions(member, send_messages=True)
                    except discord.Forbidden:
                        logger.warning(
                            f"Unable to untmute {member.name} ({member.id}) in {channel.name} ({channel.id})."
                        )
        if voice_mutes:
            for voice_mute in voice_mutes:
                channel = guild.get_channel(voice_mute.channel_snowflake)
                if channel and member.voice and member.voice.mute:
                    await member.edit(mute=False)
        await Ban.delete(**kwargs)
        await Flag.delete(**kwargs)
        await TextMute.delete(**kwargs)
        await VoiceMute.delete(**kwargs)

    @classmethod
    def add_invincible_member(cls, guild_snowflake: int, member_snowflake: int):
        cls.invincible_members[(guild_snowflake, member_snowflake)] = True

    @classmethod
    def get_invincible_members(cls):
        return cls.invincible_members

    @classmethod
    def remove_invincible_member(cls, guild_snowflake: int, member_snowflake: int):
        cls.invincible_members.pop((guild_snowflake, member_snowflake), None)

    @classmethod
    def toggle_enabled(cls):
        cls.state = not cls.state
        return cls.state
