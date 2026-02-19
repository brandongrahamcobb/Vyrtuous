"""!/bin/python3
message_service.py  The purpose of this program is to handle messages in Discord.

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

import asyncio
from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.logger import logger


class MessageService:
    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def send_message(
        self,
        source: Union[commands.Context, discord.Interaction, discord.Message],
        *,
        content: str,
        file: discord.File,
        embed: discord.Embed,
        allowed_mentions: discord.AllowedMentions = discord.AllowedMentions.none(),
        ephemeral: bool = True,
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
        allowed_mentions=None,
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


class PaginatorService:
    NAV_EMOJIS = {"\u2b05\ufe0f": -1, "\u27a1\ufe0f": 1}

    def __init__(self, bot):
        self.__bot = bot
        self.current_page = 0
        self.message = None
        self._reaction_lock = asyncio.Lock()

    async def start(self, channel, pages, *, timeout=60):
        embed = self.get_current_embed(channel=channel, pages=pages)
        if isinstance(channel, discord.Interaction):
            if not channel.response.is_done():
                await channel.response.send_message(embed=embed)
            self.message = await channel.original_response()
        elif isinstance(channel, discord.Message):
            self.message = await channel.reply(embed=embed)
        else:
            self.message = await channel.send(embed=embed)
        for emoji in self.NAV_EMOJIS:
            await self.message.add_reaction(emoji)
        self.__bot.loop.create_task(
            self.wait_for_reactions(channel=channel, pages=pages, timeout=timeout)
        )
        return self.message

    def get_current_embed(self, channel, pages):
        embed = pages[self.current_page].copy()
        total_pages = len(pages)
        label = "page"
        embed.set_footer(
            text=f"{label} {self.current_page + 1}/{total_pages} • {channel.guild.name}"
        )
        return embed

    async def wait_for_reactions(self, channel, pages, timeout):
        def look(reaction, user):
            return (
                reaction.message.id == self.message.id
                and str(reaction.emoji) in self.NAV_EMOJIS
                and not user.bot
            )

        while True:
            try:
                reaction, user = await self.__bot.wait_for(
                    "reaction_add", timeout=timeout, check=look
                )
            except asyncio.TimeoutError:
                try:
                    await self.message.clear_reactions()
                except Exception as e:
                    logger.warning(str(e).capitalize())
                break
            await self.handle_reaction(channel=channel, pages=pages, reaction=reaction)
            try:
                await self.message.remove_reaction(reaction.emoji, user)
            except Exception as e:
                logger.warning(str(e).capitalize())

    async def handle_reaction(self, channel, pages, reaction):
        async with self._reaction_lock:
            action = self.NAV_EMOJIS[str(reaction.emoji)]
            if isinstance(action, int):
                self.current_page = max(
                    0, min(self.current_page + action, len(self.pages) - 1)
                )
                await self.message.edit(
                    embed=self.get_current_embed(channel=channel, pages=pages)
                )
