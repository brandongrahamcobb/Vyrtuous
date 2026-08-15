"""!/bin/python3

discord_bot.py This is essentially a stripped version of Rapptz advanced_startup.py.

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

import logging
from typing import Self, cast

import asyncpg
import discord
from discord.ext import commands

from vyrtuous.cache.registry import (
    ChannelState,
    MemberState,
    MessageHistoryState,
    PermissionState,
    Registry,
    SystemResourcesState,
    VideoChannelState,
)


class TargetIsBot(commands.CheckFailure):
    def __init__(
        self,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        self._source = ctx or interaction or message
        if self._source is None or self._source.guild is None:
            return
        super().__init__(
            message=f"You cannot execute actions on {self._source.guild.me.mention}."
        )


class DiscordBot(commands.Bot):
    _instance = None

    def __init__(
        self,
        *,
        config,
        db_pool: asyncpg.Pool,
        initial_extensions: list[str],
        logger: logging.Logger,
        **kwargs,
    ):
        DiscordBot._instance = self
        intents = discord.Intents.all()
        # intents.message_content = False
        # intents.presences = False
        super().__init__(
            command_prefix=config["discord_command_prefix"],
            help_command=None,
            intents=intents,
            **kwargs,
        )
        self.config = config
        self.db_pool = db_pool
        self.__initial_extensions = initial_extensions
        self.logger = logger
        self.registry = Registry()
        self.testing_guild_snowflake = self.config["discord_testing_guild_snowflake"]

    async def setup_hook(self) -> None:
        for ext in self.__initial_extensions:
            await self.load_extension(ext)

    def register(self) -> None:
        self.registry.register(
            (
                ChannelState(),
                MemberState(),
                MessageHistoryState(),
                PermissionState(),
                SystemResourcesState(),
                VideoChannelState(),
            )
        )

    @classmethod
    def get_instance(cls) -> Self:
        if cls._instance is None:
            raise RuntimeError("DiscordBot instance has not been created yet")
        return cast(Self, cls._instance)
