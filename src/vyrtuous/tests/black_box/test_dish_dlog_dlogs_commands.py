"""test_trole_xtrole_commands.py The purpose of this program is to black box test the team role commands.

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
"""

from typing import Optional

import pytest
import discord

from vyrtuous.inc.helpers import PRIVILEGED_AUTHOR_ID
from vyrtuous.database.roles.developer import Developer
from vyrtuous.tests.black_box.test_suite import (
    extract_embed_text,
    prepared_command_handling,
    RESET,
    YELLOW,
    RED,
    GREEN,
)
from vyrtuous.utils.emojis import EMOJIS


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "command",
    [
        (f"dish 5cb608e7-7b95-4c22-9bd4-aac414562b10 {PRIVILEGED_AUTHOR_ID}"),
        (f"dish 5cb608e7-7b95-4c22-9bd4-aac414562b10 {PRIVILEGED_AUTHOR_ID}"),
        ("dlogs unresolved 4588f2ca-1a03-4b6d-beb9-b88cbc8b2b58"),
        ("dlogs resolved 5cb608e7-7b95-4c22-9bd4-aac414562b10"),
        ("dlog 4588f2ca-1a03-4b6d-beb9-b88cbc8b2b58 overwrite Test"),
        ("dlog 4588f2ca-1a03-4b6d-beb9-b88cbc8b2b58 append test2"),
        ("dlog 4588f2ca-1a03-4b6d-beb9-b88cbc8b2b58 resolve"),
    ],
)
async def test_dish_dlog_dlogs_commands(
    bot,
    voice_channel_one,
    guild,
    privileged_author,
    prefix: Optional[str],
    command: Optional[str]
):
    developer = Developer(
        guild_snowflake=guild.id, member_snowflake=privileged_author.id
    )
    await developer.create()
    try:
        voice_channel_one.messages.clear()
        formatted = f"{command}"
        captured = await prepared_command_handling(
            author=privileged_author,
            bot=bot,
            channel=voice_channel_one,
            cog="DevCommands",
            content=formatted,
            guild=guild,
            isinstance_patch="vyrtuous.cogs.dev_commands.isinstance",
            prefix=prefix,
        )
        message = captured["message"]
        message_type = captured["type"]
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
            print(f"{GREEN}Success:{RESET} {content}")
            assert any(emoji in content for emoji in EMOJIS)
    finally:
        await developer.delete()
