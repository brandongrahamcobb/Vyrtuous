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

    async def resolve_member(self, ctx_interaction_or_message, scope):
        guild = ctx_interaction_or_message.guild
        if isinstance(scope, discord.Member):
            return scope
        if isinstance(scope, int):
            member_id = scope
        elif isinstance(scope, str):
            if scope.isdigit():
                member_id = int(scope)
            elif scope.startswith('<@') and scope.endswith('>'):
                member_id = int(scope[2:-1].replace('!', ''))
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
            raise commands.BadArgument(f'Member `{scope}` not found in {ctx_interaction_or_message.guild.name}.')

