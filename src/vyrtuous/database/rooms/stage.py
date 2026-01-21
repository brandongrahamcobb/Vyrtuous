"""stage.py The purpose of this program is to provide the Stage utility class.

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

from datetime import datetime
from typing import Optional
import time

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.rooms.room import Room
from vyrtuous.utils.emojis import get_random_emoji


class Stage(Room):

    ACT = "stage"
    PLURAL = "Stages"
    SINGULAR = "Stage"
    UNDO = "stage"
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "expires_in",
        "guild_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "expired",
        "reason",
        "target",
        "updated_at",
    ]
    TABLE_NAME = "active_stages"

    def __init__(
        self,
        channel_snowflake: int,
        expires_in: Optional[datetime],
        guild_snowflake: int,
        reason: Optional[str] = "No reason provided.",
        target: Optional[str] = "room",
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.expires_in = expires_in
        self.guild_snowflake = guild_snowflake
        self.reason = reason
        self.target = target
        self.updated_at = updated_at

    async def send_stage_ask_to_speak_message(
        self, join_log: dict[int, discord.Member], member: discord.Member
    ):
        bot = DiscordBot.get_instance()
        now = time.time()
        join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
        if len(join_log[member.id]) < 1:
            join_log[member.id].append(now)
            embed = discord.Embed(
                title=f"{get_random_emoji()} {self.channel_snowflake} — Stage Mode",
                description=f"Ends <t:{int(self.expires_in.timestamp())}:R>",
                color=discord.Color.green(),
            )
            embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
            await bot.get_channel(self.channel_snowflake).send(embed=embed)
