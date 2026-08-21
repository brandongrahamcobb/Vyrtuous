"""!/bin/python3
mock_discord_bot.py The purpose of this program is to support integration testing for Vyrtuous.

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

from unittest.mock import AsyncMock, Mock

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
from vyrtuous.inc.helpers import DISCORD_COGS, PATH_LOG
from vyrtuous.system.config import Config
from vyrtuous.system.logger import logger, setup_logging


class MockBot(commands.Bot):
    def __init__(self, config: Config, db_pool: asyncpg.Pool):
        intents = discord.Intents.all()
        intents.message_content = False
        setup_logging(config, PATH_LOG)
        self.config = config
        self.db_pool = db_pool
        self._guilds = {}
        self.registry = Registry()
        self._tree = AsyncMock()
        self._tree.sync = AsyncMock()
        self._tree.add_command = Mock()
        self._tree.remove_command = AsyncMock()
        self._tree.copy_global_to = AsyncMock()
        self.logger = logger
        super().__init__(
            command_prefix=config.get("discord_command_prefix", "!"),
            help_command=None,
            intents=intents,
        )

    def register(self):
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
    def get_instance(cls):
        return cls()

    def get_guild(self, guild_snowflake: int):
        return self._guilds.get(guild_snowflake, None)

    def get_channel(self, channel_snowflake: int):
        for guild in self._guilds.values():
            for channel in guild.channels:
                if channel.id == channel_snowflake:
                    return channel
        return None

    def get_user(self, member_snowflake: int):
        for guild_snowflake, guild in self._guilds.items():
            for member in guild.members:
                if member.id == member_snowflake:
                    return member

    async def setup_hook(self):
        for cog in DISCORD_COGS:
            if cog != "vyrtuous.utils.scheduled_tasks":
                await self.load_extension(cog)

    @property
    def tree(self):
        return self._tree

    @property
    def guilds(self):
        return list(self._guilds.values())
