''' member_service.py The purpose of this program is to provide the member_service module.
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
from discord.ext import commands
from typing import Optional, Union
from vyrtuous.utils.setup_logging import logger

import discord

class MemberService:

    async def resolve_member(self, ctx_interaction_or_message, member_str):
        guild = ctx_interaction_or_message.guild
        if isinstance(member_str, discord.Member):
            return member_str
        if isinstance(member_str, int):
            member_id = member_str
        elif isinstance(member_str, str):
            if member_str.isdigit():
                member_id = int(member_str)
            elif member_str.startswith('<@') and member_str.endswith('>'):
                member_id = int(member_str[2:-1].replace('!', ''))
            else:
                raise commands.BadArgument('Invalid member identifier')
        else:
            raise commands.BadArgument('Invalid member identifier')
        member = guild.get_member(member_id)
        if member:
            return member
        try:
            member = await guild.fetch_member(member_id)
            return member
        except discord.NotFound:
            raise commands.BadArgument(f'Member `{member_str}` not found in {ctx_interaction_or_message.guild.name}.')

