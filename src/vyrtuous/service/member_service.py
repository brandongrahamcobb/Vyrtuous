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
from typing import Optional, Union
from vyrtuous.utils.setup_logging import logger

import discord

class MemberService:

    async def resolve_member(self, ctx_interaction_or_message, value: Optional[Union[int, str, discord.Member]]) -> Optional[discord.Member]:
        try:
            if isinstance(value, discord.Member):
                logger.debug(f"Direct member: {value.id}")
                return value
            if isinstance(value, int):
                m = ctx_interaction_or_message.guild.get_member(value)
                if not m:
                    try: m = await ctx_interaction_or_message.guild.fetch_member(value)
                    except discord.NotFound: m = None
                if m:
                    logger.debug(f"Resolved member by int ID: {m.id}")
                    return m
            if isinstance(value, str):
                if value.isdigit():
                    mid = int(value)
                    m = ctx_interaction_or_message.guild.get_member(mid)
                    if not m:
                        try: m = await ctx_interaction_or_message.guild.fetch_member(mid)
                        except discord.NotFound: m = None
                    if m:
                        logger.debug(f"Resolved member by str ID: {m.id}")
                        return m
                if value.startswith('<@') and value.endswith('>'):
                    mid = int(value[2:-1].replace('!', ''))
                    m = ctx_interaction_or_message.guild.get_member(mid)
                    if not m:
                        try: m = await ctx_interaction_or_message.guild.fetch_member(mid)
                        except discord.NotFound: m = None
                    if m:
                        logger.debug(f"Member mention resolved: {m.id}")
                        return m
                    return m
        except Exception as e:
            logger.warning(f"Member resolution error: {e}")
        return None
