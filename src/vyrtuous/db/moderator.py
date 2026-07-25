"""!/bin/python3
moderator.py The purpose of this program is to provide the moderator database model.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from dataclasses import dataclass, field
from datetime import datetime, timezone

from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot


@dataclass(frozen=True)
class Moderator:
    __tablename__ = "moderators"
    identifier = "mod"
    channel_snowflake: int
    guild_snowflake: int
    member_snowflake: int
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))


class NotModerator(commands.CheckFailure):
    def __init__(
        self,
        channel_snowflake: int | None,
        guild_snowflake: int,
    ):
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        if channel_snowflake is None:
            message: str = (
                f"You lack sufficient permissions of a moderator in the requested server ({guild.name})."
            )
            super().__init__(message)
        else:
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                raise commands.ChannelNotFound(str(channel_snowflake))
            message: str = (
                f"You lack sufficient permissions of a moderator in the requested channel ({channel.mention})."
            )
            super().__init__(message)


class NotAppModerator(app_commands.CheckFailure):
    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
    ):
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        channel = guild.get_channel(channel_snowflake)
        if channel is None:
            raise commands.ChannelNotFound(str(channel_snowflake))
        message: str = (
            f"You lack sufficient permissions of a moderator in the requested channel ({channel.mention})."
        )
        super().__init__(message)
