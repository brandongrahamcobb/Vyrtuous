from typing import Optional, Union

import discord
from vyrtuous.utils.setup_logging import logger

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
