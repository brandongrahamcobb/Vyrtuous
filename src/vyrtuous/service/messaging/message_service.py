"""discord_message_service.py  The purpose of this program is to handle messages in Discord.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from discord.ext import commands
from typing import Union
import discord

from vyrtuous.bot.discord_bot import DiscordBot


class MessageService:

    def __init__(self, bot: DiscordBot, db_pool):
        self.bot = bot

    async def send_message(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
        *,
        content: str = None,
        file: discord.File = None,
        embed: discord.Embed = None,
        allowed_mentions: discord.AllowedMentions = discord.AllowedMentions.none(),
        ephemeral: bool = None
    ):
        if isinstance(source, commands.Context):
            can_send = (
                source.guild
                and isinstance(source.channel, discord.abc.GuildChannel)
                and source.channel.permissions_for(source.guild.me).send_messages
            )
            if can_send:
                try:
                    return await self._send_message(
                        lambda **kw: source.reply(**kw),
                        content=content,
                        file=file,
                        embed=embed,
                        allowed_mentions=allowed_mentions,
                    )
                except discord.HTTPException as e:
                    if getattr(e, "code", None) == 50035:
                        return await self._send_message(
                            lambda **kw: source.send(**kw),
                            content=content,
                            file=file,
                            embed=embed,
                            allowed_mentions=allowed_mentions,
                        )
                    else:
                        raise
            else:
                return await self.send_dm(
                    source.author,
                    content=content,
                    file=file,
                    embed=embed,
                    allowed_mentions=allowed_mentions,
                )
        elif isinstance(source, discord.Interaction):
            if ephemeral is None:
                ephemeral = True
            if source.response.is_done():
                return await self._send_message(
                    lambda **kw: source.followup.send(**kw, ephemeral=ephemeral),
                    content=content,
                    file=file,
                    embed=embed,
                    allowed_mentions=allowed_mentions,
                )
            else:
                await self._send_message(
                    lambda **kw: source.response.send_message(
                        **kw, ephemeral=ephemeral
                    ),
                    content=content,
                    file=file,
                    embed=embed,
                    allowed_mentions=allowed_mentions,
                )
                return await source.original_response()
        elif isinstance(source, discord.Message):
            return await self._send_message(
                lambda **kw: source.reply(**kw),
                content=content,
                file=file,
                embed=embed,
                allowed_mentions=allowed_mentions,
            )
        else:
            raise TypeError(
                "Expected commands.Context, discord.Interaction, or discord.Message"
            )

    async def _send_message(
        self, send_func, *, content=None, file=None, embed=None, allowed_mentions=None
    ):
        kwargs = {}
        if content:
            kwargs["content"] = content
        if file:
            kwargs["file"] = file
        if embed:
            kwargs["embed"] = embed
        if allowed_mentions is not None:
            kwargs["allowed_mentions"] = allowed_mentions
        return await send_func(**kwargs)

    async def send_dm(
        self,
        user: discord.abc.User,
        *,
        content=None,
        file=None,
        embed=None,
        allowed_mentions=None
    ):
        dm_channel = user.dm_channel
        if dm_channel is None:
            dm_channel = await user.create_dm()
        return await self._send_message(
            lambda **kw: dm_channel.send(**kw),
            content=content,
            file=file,
            embed=embed,
            allowed_mentions=allowed_mentions,
        )
