''' test_log_logs_mlog_commands.py The purpose of this program is to black box test the logging related commands.

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
from vyrtuous.tests.black_box.test_suite import *
from vyrtuous.enhanced_members.administrator import Administrator
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,action,target_type,text_ref, channel_ref,member_ref",
    [
        ("mlog", "create", "general", True, False, False),
        ("log", None, None, True, False, False),
        ("log", None, None, True, False, False),
        ("logs", "all", None, False, False, False),
        ("logs", None, None, False, False, False),
        ("mlog", "modify", "general", True, False, False),
        ("mlog", "delete", "general", True, False, False),
        ("mlog", "create", "channel", True, True, False),
        ("mlog", "modify", "channel", True, True, False),
        ("mlog", "delete", None, True, False, False),
        ("mlog", "create", "member", True, False, True),
        ("mlog", "modify", "member", True, False, True),
        ("mlog", "delete", None, True, False, False)
    ]
)

async def test_log_logs_mlog_command(bot, text_channel, voice_channel_one, guild, not_privileged_author, privileged_author, prefix: Optional[str], role, command: Optional[str], action: Optional[str], target_type: Optional[str], channel_ref, member_ref, text_ref):
    administrator = Administrator(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflakes=[role.id])
    await administrator.grant()
    try:
        text_channel.messages.clear() 
        text_channel_value = text_channel.mention
        voice_channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
        privileged_author_value = privileged_author.mention if member_ref else privileged_author.name
        not_privileged_author_value = not_privileged_author.mention if member_ref else not_privileged_author.name
        if command == "mlog":
            match target_type:
                case "general":
                    formatted = f"{command} {text_channel.id} {action} {target_type}".strip()
                case "channel":
                    formatted = f"{command} {text_channel.id} {action} {target_type} {voice_channel_one.id}".strip()
                case "member":
                    formatted = f"{command} {text_channel.id} {action} {target_type} {not_privileged_author.id}".strip()
                case None:
                    formatted = f"{command} {text_channel.id} {action}".strip()
        if command == "log":
            formatted = f"{command} {text_channel.id}".strip()
        if command == "logs" and action:
            formatted = f"{command} {action}".strip()
        elif command == "logs":
            formatted = f"{command}".strip()
        captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=text_channel, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
        message = captured['message']
        message_type = captured['type']
        if isinstance(message, discord.Embed):
            content = extract_embed_text(message)
        elif isinstance(message, discord.File):
            content = message.filename
        else:
            content = message
        if message_type == "error":
            print(f"{RED}Error:{RESET} {content}")
        if message_type == "warning":
            print(f"{YELLOW}Warning:{RESET} {content}")
        if message_type == "success":
            # print(f"{GREEN}Success:{RESET} {content}")
            # assert any(emoji in content for emoji in Emojis.EMOJIS)
            if text_ref:
                assert text_channel_value in content
            # match target_type:
            #     case "channel":
            #         assert voice_channel_value in content
            #     case "member":
            #         assert not_privileged_author_value in content
            # if channel_ref:
            #     assert text_channel_value in content
    finally:
        await administrator.revoke()