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

import discord

from vyrtuous.bot.discord_bot import DiscordBot


class Paginator:
    NAV_EMOJIS = {"\u2b05\ufe0f": -1, "\u27a1\ufe0f": 1}

    def __init__(self):
        self.current_page = 0
        self._reaction_lock = asyncio.Lock()

    async def start(
        self, channel, pages, *, ephemeral=False, timeout=60
    ) -> discord.Message:
        bot: DiscordBot = DiscordBot.get_instance()
        embed = self.get_current_embed(channel=channel, pages=pages)
        if isinstance(channel, discord.Interaction):
            if not channel.response.is_done():
                await channel.response.send_message(embed=embed, ephemeral=ephemeral)
            message = await channel.original_response()
        elif isinstance(channel, discord.Message):
            message = await channel.reply(embed=embed)
        else:
            message = await channel.send(embed=embed)
        for emoji in self.NAV_EMOJIS:
            await message.add_reaction(emoji)
        bot.loop.create_task(
            self.wait_for_reactions(
                channel=channel, message=message, pages=pages, timeout=timeout
            )
        )
        return message

    def get_current_embed(self, channel, pages) -> discord.Embed:
        embed = pages[self.current_page].copy()
        total_pages = len(pages)
        label = "page"
        embed.set_footer(
            text=f"{label} {self.current_page + 1}/{total_pages} • {channel.guild.name}"
        )
        return embed

    async def wait_for_reactions(self, channel, message, pages, timeout) -> None:
        bot: DiscordBot = DiscordBot.get_instance()

        def look(reaction, user):
            return (
                reaction.message.id == message.id
                and str(reaction.emoji) in self.NAV_EMOJIS
                and not user.bot
            )

        while True:
            try:
                reaction, user = await bot.wait_for(
                    "reaction_add", timeout=timeout, check=look
                )
            except asyncio.TimeoutError:
                try:
                    await message.clear_reactions()
                except Exception as e:
                    bot.logger.warning(str(e).capitalize())
                break
            await self.handle_reaction(
                channel=channel, message=message, pages=pages, reaction=reaction
            )
            try:
                await message.remove_reaction(reaction.emoji, user)
            except Exception as e:
                bot.logger.warning(str(e).capitalize())

    async def handle_reaction(self, channel, message, pages, reaction) -> None:
        async with self._reaction_lock:
            action = self.NAV_EMOJIS[str(reaction.emoji)]
            if isinstance(action, int):
                self.current_page = max(
                    0, min(self.current_page + action, len(pages) - 1)
                )
                await message.edit(
                    embed=self.get_current_embed(channel=channel, pages=pages)
                )
