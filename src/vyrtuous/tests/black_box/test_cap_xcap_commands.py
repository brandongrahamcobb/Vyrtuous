# ''' test_cap_xcap_commands.py The purpose of this program is to black box test the Cap commands.
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
# import pytest

# def generate_cap_test_cases():
#     durations = ['0', '1', '24', '48']
#     moderation_types = ['ban', 'vmute', 'tmute']
#     cases = []
#     for mod_type in moderation_types:
#         for duration in durations:
#             cases.append((f"cap {{voice_channel_one_id}} {mod_type} {duration}", mod_type, True))
#             cases.append((f"xcap {{voice_channel_one_id}} {mod_type}", mod_type, True))
#     return cases

# @pytest.mark.asyncio
# @pytest.mark.parametrize("command,moderation_type,channel_ref", generate_cap_test_cases())
# async def test_cap_commands(bot, voice_channel_one, guild, privileged_author, prefix: Optional[str], role, command: Optional[str], moderation_type, channel_ref):
#     administrator = Administrator(guild_snowflake=guild.id, member_snowflake=privileged_author.id, role_snowflakes=[role.id])
#     await administrator.grant()
#     try:
#         channel_value = voice_channel_one.mention if channel_ref else voice_channel_one.name
#         match moderation_type:
#             case "tmute":
#                 moderation_type = "text_mute"
#             case "vmute":
#                 moderation_type = "voice_mute"
#             case "untmute":
#                 moderation_type = "untext_mute"
#             case "unvmute":
#                 moderation_type = "unvoice_mute"
#         voice_channel_one.messages.clear()
#         formatted = command.format(
#             voice_channel_one_id=voice_channel_one.id
#         )
#         captured = await prepared_command_handling(author=privileged_author, bot=bot, channel=voice_channel_one, cog="AdminCommands", content=formatted, guild=guild, isinstance_patch="vyrtuous.cogs.admin_commands.isinstance", prefix=prefix)
#         message = captured['message']
#         message_type = captured['type']
#         if isinstance(message, discord.Embed):
#             content = extract_embed_text(message)
#         elif isinstance(message, discord.File):
#             content = message.filename
#         else:
#             content = message
#         if message_type == "error":
#             print(f"{RED}Error:{RESET} {content}")
#         if message_type == "warning":
#             print(f"{YELLOW}Warning:{RESET} {content}")
#         if message_type == "success":
#             assert any(emoji in content for emoji in EMOJIS)
#             assert moderation_type in content
#             if channel_ref:
#                 assert channel_value in content
#     finally:
#         await administrator.revoke()
