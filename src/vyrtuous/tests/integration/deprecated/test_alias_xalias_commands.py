# ''' test_alias_xalias_commands.py The purpose of this program is to black box test the Aliases module.
#     Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
# '''
# from typing import Optional
# from vyrtuous.inc.helpers import *
# from vyrtuous.tests.black_box.test_suite import *
# from vyrtuous.database.roles.administrator import Administrator
# from vyrtuous.utils.emojis import get_random_emoji, EMOJIS
# import discord
# import pytest

# @pytest.mark.asyncio
# @pytest.mark.parametrize(
#     "command,alias_type,alias_name,channel_ref,role_ref",
#     [
#         ("alias", "ban", "testban", True, False),
#         ("xalias", None, "testban", True, False),
#         ("alias", "unban", "testunban", True, False),
#         ("xalias", None, "testunban", True, False),
#         ("alias", "voice_mute", "testmute", True, False),
#         ("xalias", None, "testmute", True, False),
#         ("alias", "unvoice_mute", "testunmute", True, False),
#         ("xalias", None, "testunmute", True, False),
#         ("alias", "flag", "testflag", True, False),
#         ("xalias", None, "testflag",  True, False),
#         ("alias", "unflag", "testunflag", True, False),
#         ("xalias", None, "testunflag", True, False),
#         ("alias", "vegan", "testvegan", True, False),
#         ("xalias", None, "testvegan", True, False),
#         ("alias", "carnist", "testcarnist", True, False),
#         ("xalias", None, "testcarnist", True, False),
#         ("alias", "tmute", "testtmute", True, False),
#         ("xalias", None, "testtmute", True, False),
#         ("alias", "untmute", "testuntmute", True, False),
#         ("xalias", None, "testuntmute", True, False),
#         ("alias", "role", "testrole", True, True),
#         ("xalias", "role", "testrole", True, True),
#         ("alias", "unrole", "testunrole", True, True),
#         ("xalias", "unrole", "testunrole", True, True)
#     ]
# )

# async def test_alias_xalias_command(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], command: Optional[str], role, alias_type, alias_name, channel_ref, role_ref):
#     administrator = Administrator(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflakes=[role.id])
#     await administrator.create()
#     try:
#         channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
#         role_value = role.mention if role_ref else role.name
#         voice_channel_one.messages.clear()
#         if channel_ref and role_ref and command == "alias":
#             formatted = f"{command} {alias_type} {alias_name} {voice_channel_one.id} {ROLE_ID}"
#         elif channel_ref and command == "alias":
#             formatted = f"{command} {alias_type} {alias_name} {voice_channel_one.id}"
#         else:
#             formatted = f"{command} {alias_name}"
#         captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
#         message = captured['message']
#         message_type = captured['type']
#         if isinstance(message, discord.Embed):
#             content = extract_embed_text(message)
#         else:
#             content = message
#         if message_type == "error":
#             print(f"{RED}Error:{RESET} {content}")
#         if message_type == "warning":
#             print(f"{YELLOW}Warning:{RESET} {content}")
#         if message_type == "success":
#             # print(f"{GREEN}Success:{RESET} {content}")
#             if alias_type:
#                 assert alias_type in content and alias_name in content and channel_value in content
#             if role_ref:
#                 assert role_value in content
#     finally:
#         await administrator.delete()
