''' test_mlog_command.py The purpose of this program is to black box test the mlog command.
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
from vyrtuous.tests.test_admin_helpers import admin_cleanup, admin_initiation, prepared_command_handling
from vyrtuous.tests.test_suite import *
from vyrtuous.utils.emojis import Emojis
import discord
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,action,target_type,target_id,channel_ref,member_ref",
    [
        ("mlog", "create", "general", None, None, None),
        ("mlog", "modify", "general", None, None, None),
        ("mlog", "delete", "general", None, None, None),
        ("mlog", "create", "channel", "voice_channel_one_id", True, None),
        ("mlog", "modify", "channel", "voice_channel_one_id", True, None),
        ("mlog", "delete", None, "voice_channel_one_id", True, None),
        ("mlog", "create", "member", "voice_channel_one_id", None, None),
        ("mlog", "modify", "member", "voice_channel_one_id", None, None),
        ("mlog", "delete", None, "voice_channel_one_id", None, None)
    ]
)

async def test_mlog_command(bot, text_channel, voice_channel_one, guild, privileged_author, not_privileged_author, prefix: Optional[str], command: Optional[str], action: Optional[str], target_type: Optional[str], target_id: Optional[str], channel_ref, member_ref):
    await admin_initiation(guild.id, privileged_author.id)
    try:
        text_channel.messages.clear() 
        target_obj = {
            "voice_channel_one_id": voice_channel_one.id
        }.get(target_id, None)
        if target_obj:
            target_str = target_obj
            formatted = f"{command} {text_channel.id} {action} {target_type} {target_str}".strip()
        else:
            formatted = f"{command} {text_channel.id} {action} {target_type}".strip()
        await prepared_command_handling(author=privileged_author, bot=bot, channel=text_channel, content=formatted, guild=guild, prefix=prefix)
        response = text_channel.messages[0]
        print(response)
        channel_value = text_channel.mention if channel_ref else text_channel.name
        self_member_value = privileged_author.mention if member_ref else privileged_author.name
        dummy_member_value = not_privileged_author.mention if member_ref else not_privileged_author.name
        assert any(emoji in response["content"] for emoji in Emojis.EMOJIS)
        if channel_ref:
            assert any(val in response["content"] for val in [channel_value])
        if member_ref == "not_privileged_author":
            assert any(val in response["content"] for val in [dummy_member_value])
        if member_ref == "privileged_author":
            assert any(val in response["content"] for val in [self_member_value])
        text_channel.messages.clear() 
    finally:
        await admin_cleanup(guild.id, privileged_author.id)