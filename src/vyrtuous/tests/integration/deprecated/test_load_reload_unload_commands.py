# """test_load_reload_unload_commands.py The purpose of this program is to black box test the module management commands.

# Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# """

# from typing import Optional

# import pytest

# from vyrtuous.inc.helpers import DISCORD_COGS
# from vyrtuous.tests.black_box.test_suite import (
#     extract_embed_text,
#     prepared_command_handling,
#     RESET,
#     YELLOW,
#     RED,
#     GREEN,
# )
# from vyrtuous.utils.emojis import EMOJIS


# @pytest.mark.asyncio
# @pytest.mark.parametrize(
#     "permission,command,should_warn",
#     [
#         ("Developer", f"{action} {cog}", False)
#         for action in ("load", "reload", "unload")
#         for cog in DISCORD_COGS
#         if cog != "vyrtuous.cogs.dev_commands"
#     ],
#     indirect=["permission"],
# )
# async def test_load_reload_unload_command(
#     bot,
#     command: Optional[str],
#     permission,
#     prefix: Optional[str],
#     privileged_author,
#     text_channel,
#     should_warn,
#     guild,
# ):
#     formatted = command
#     captured = await prepared_command_handling(
#         author=privileged_author,
#         bot=bot,
#         channel=text_channel,
#         content=formatted,
#         guild=guild,
#         highest_role=permission,
#         prefix=prefix,
#     )
#     message = captured[0]["message"]
#     message_type = captured[0]["type"]
#     if message.embeds:
#         embed = message.embeds[0]
#         content = extract_embed_text(embed)
#     elif message.embed:
#         content = extract_embed_text(message.embed)
#     else:
#         content = message.content
#     if message_type == "error":
#         print(f"{RED}Error:{RESET} {content}")
#     if message_type == "warning":
#         print(f"{YELLOW}Warning:{RESET} {content}")
#         if should_warn:
#             assert True
#     if message_type == "success":
#         print(f"{GREEN}Success:{RESET} {content}")
#         assert any(emoji in content for emoji in EMOJIS)
