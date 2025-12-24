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
from typing import Optional
from vyrtuous.tests.black_box.test_suite import bot, config, guild, not_privileged_author, prepared_command_handling, prefix, privileged_author, role, text_channel, voice_channel_one
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,action,target_type,target_id,channel_ref,member_ref",
    [
        ("mlog", "create", "general", None, None, None),
        ("log", None, None, None, True, False),
        ("log", None, None, None, True, False),
        ("logs", "all", None, None, False, False),
        ("logs", None, None, None, True, False),
        ("mlog", "modify", "general", None, None, None),
        ("mlog", "delete", "general", None, None, None),
        ("mlog", "create", "channel", True, True, None),
        ("mlog", "modify", "channel", True, True, None),
        ("mlog", "delete", None, True, True, None),
        ("mlog", "create", "member", True, None, None),
        ("mlog", "modify", "member", True, None, None),
        ("mlog", "delete", None, True, None, None)
    ]
)

async def test_log_logs_mlog_command(bot, text_channel, voice_channel_one, guild, privileged_author, prefix: Optional[str], role, command: Optional[str], action: Optional[str], target_type: Optional[str], target_id: Optional[str], channel_ref, member_ref):
    await Administrator.grant(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflake=role.id)
    try:
        text_channel.messages.clear() 
        if command == "mlog":
            if target_id:
                formatted = f"{command} {text_channel.id} {action} {target_type} {voice_channel_one.id}".strip()
            else:
                formatted = f"{command} {text_channel.id} {action} {target_type}".strip()
        if command in ("log", "logs"):
            if channel_ref:
                formatted = f"{command} {text_channel.id}".strip()
            else:
                formatted = f"{command} {action}".strip()
        await prepared_command_handling(author=privileged_author, bot=bot, channel=text_channel, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
        response = text_channel.messages[0]
        channel_value = text_channel.mention if channel_ref else text_channel.name
        privileged_author_value = privileged_author.mention if member_ref else privileged_author.name
        not_privileged_author_value = not_privileged_author.mention if member_ref else not_privileged_author.name
        if command in ("log", "mlog"):
            assert any(emoji in response["content"] for emoji in Emojis.EMOJIS)
            if channel_ref:
                assert any(val in response["content"] for val in [channel_value])
        else:
            assert any(emoji in response["embed"].title for emoji in Emojis.EMOJIS)
        if member_ref:
            assert any(val in response["content"] for val in [privileged_author_value])
    finally:
        await Administrator.revoke(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflake=role.id)