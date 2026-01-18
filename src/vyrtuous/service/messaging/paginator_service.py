"""paginator.py  The purpose of this program is to handle multiple page embeds in Discord.

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
"""

import asyncio

import discord

from vyrtuous.service.logging_service import logger


class Paginator:

    NAV_EMOJIS = {"\u2b05\ufe0f": -1, "\u27a1\ufe0f": 1}

    def __init__(self, bot, channel_source, pages, state=None, *, timeout=60):
        self.bot = bot
        self.channel_source = channel_source
        self.pages = pages
        self.state = state
        self.current_page = 0
        self.timeout = timeout
        self.message = None
        self._reaction_lock = asyncio.Lock()

    async def start(self):
        embed = self.get_current_embed()
        if isinstance(self.channel_source, discord.Interaction):
            if not self.channel_source.response.is_done():
                await self.channel_source.response.send_message(embed=embed)
            self.message = await self.channel_source.original_response()
        elif isinstance(self.channel_source, discord.Message):
            self.message = await self.channel_source.reply(embed=embed)
        else:
            self.message = await self.channel_source.send(embed=embed)
        for emoji in self.NAV_EMOJIS:
            await self.message.add_reaction(emoji)
        self.bot.loop.create_task(self.wait_for_reactions())
        return self.message

    def get_current_embed(self):
        embed = self.pages[self.current_page].copy()
        total_pages = len(self.pages)
        label = "page"
        embed.set_footer(
            text=f"{label} {self.current_page + 1}/{total_pages} • {self.channel_source.guild.name}"
        )
        return embed

    async def wait_for_reactions(self):
        def look(reaction, user):
            return (
                reaction.message.id == self.message.id
                and str(reaction.emoji) in self.NAV_EMOJIS
                and not user.bot
            )

        while True:
            try:
                reaction, user = await self.bot.wait_for(
                    "reaction_add", timeout=self.timeout, check=look
                )
            except asyncio.TimeoutError:
                try:
                    await self.message.clear_reactions()
                except Exception as e:
                    logger.warning(str(e).capitalize())
                break
            await self.handle_reaction(reaction)
            try:
                await self.message.remove_reaction(reaction.emoji, user)
            except Exception as e:
                logger.warning(str(e).capitalize())

    async def handle_reaction(self, reaction):
        async with self._reaction_lock:
            action = self.NAV_EMOJIS[str(reaction.emoji)]
            if isinstance(action, int):
                self.current_page = max(
                    0, min(self.current_page + action, len(self.pages) - 1)
                )
                await self.message.edit(embed=self.get_current_embed())
