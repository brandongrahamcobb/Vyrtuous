''' test_xrmute_command.py The purpose of this program is to black box test the unroom-mute command.

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
from vyrtuous.utils.coordinator import Coordinator
from vyrtuous.utils.emojis import Emojis
import pytest

@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command,member_ref",
    [
        ("xrmute", "True"),
        ("xrmute", "True")
    ]
)

async def test_xrmute_command(bot, voice_channel_one, guild, privileged_author, not_privileged_author, prefix: Optional[str], command: Optional[str], member_ref):
    coordinator = Coordinator(channel_snowflake=voice_channel_one.id, guild_snowflake=guild.id, member_snowflake=privileged_author.id)
    await coordinator.grant()
    try:
        voice_channel_one.messages.clear() 
        formatted = f"{command} {voice_channel_one.id}"
        captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="CoordinatorCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.coordinator_commands.isinstance", prefix=prefix)
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
            assert any(emoji in content for emoji in Emojis.EMOJIS)
    finally:
        await coordinator.revoke()