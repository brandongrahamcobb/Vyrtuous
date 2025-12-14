from typing import Optional, Union

import discord

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
