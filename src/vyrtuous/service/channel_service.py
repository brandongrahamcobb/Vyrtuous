''' channel_service.py The purpose of this program is to provide the channel_service module.
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
import discord
from typing import Optional, Union
from vyrtuous.utils.setup_logging import logger

class ChannelService:
    
    async def resolve_channel(self, ctx_interaction_or_message, value: Optional[Union[int, str, discord.TextChannel, discord.VoiceChannel]]) -> Optional[Union[discord.TextChannel, discord.VoiceChannel]]:
        try:
            if isinstance(value, (discord.TextChannel, discord.VoiceChannel)):
                logger.debug(f"Direct channel: {value.id}")
                return value
            if isinstance(value, int):
                c = ctx_interaction_or_message.guild.get_channel(value)
                if isinstance(c, (discord.TextChannel, discord.VoiceChannel)):
                    logger.debug(f"Resolved channel by int ID: {c.id}")
                    return c
            if isinstance(value, str):
                if value.isdigit():
                    cid = int(value)
                    c = ctx_interaction_or_message.guild.get_channel(cid)
                    if isinstance(c, (discord.TextChannel, discord.VoiceChannel)):
                        logger.debug(f"Resolved channel by str ID: {c.id}")
                        return c
                if value.startswith('<#') and value.endswith('>'):
                    cid = int(value[2:-1])
                    c = ctx.guild.get_channel(cid)
                    if c:
                        logger.debug(f"Channel mention resolved: {c.id}")
                        return c
        except Exception as e:
            logger.warning(f"Channel resolution error: {e}")
        return ctx_interaction_or_message.channel
