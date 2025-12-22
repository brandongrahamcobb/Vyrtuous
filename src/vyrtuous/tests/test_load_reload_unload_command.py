''' test_backup_command.py The purpose of this program is to black box test the backup command.
    Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
'''
from discord.ext.commands import view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import PropertyMock, patch
from vyrtuous.inc.helpers import *
from vyrtuous.tests.test_dev_helpers import dev_cleanup, dev_initiation
from vyrtuous.tests.test_suite import bot, config, guild, prepared_command_handling, prefix, privileged_author, voice_channel_one, voice_channel_two
from vyrtuous.utils.emojis import Emojis
import asyncio
import discord
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        f"{action} {cog}"
        for action in ("load", "reload", "unload")
        for cog in DISCORD_COGS
        if cog != "vyrtuous.cogs.dev_commands"
    ]
)

async def test_clear_command(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], command: Optional[str]):
    await dev_initiation(guild.id, privileged_author.id)
    try:
        voice_channel_one.messages.clear() 
        await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="DevCommands", content=command, guild=guild, isinstance_patch="vyrtuous.cogs.dev_commands.isinstance", prefix=prefix)
        response = voice_channel_one.messages[0]["content"]
        assert '\N{OK HAND SIGN}' in response 
    finally:
        await dev_cleanup(guild.id, privileged_author.id)