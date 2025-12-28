
''' discord_message_service.py  The purpose of this program is to handle messages in Discord.

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
import asyncio
import discord

class Paginator:

    NAV_EMOJIS = {
        "\u2b05\ufe0f": -1,
        "\u27a1\ufe0f": 1
    }

    def __init__(self, bot, ctx_or_interaction, pages, state=None, *, timeout=60):
        self.bot = bot
        self.ctx_or_interaction = ctx_or_interaction
        self.pages = pages
        self.state = state
        self.current_page = 0
        self.timeout = timeout
        self.message = None

    async def start(self):
        embed = self.pages[self.current_page]
        if isinstance(self.ctx_or_interaction, discord.Interaction):
            await self.ctx_or_interaction.response.defer()
            self.message = await self.ctx_or_interaction.followup.send(embed=embed)
        else:
            self.message = await self.ctx_or_interaction.send(embed=embed)
        for emoji in self.NAV_EMOJIS:
            await self.message.add_reaction(emoji)
        self.bot.loop.create_task(self.wait_for_reactions())
        return self.message

    async def wait_for_reactions(self):
        def check(reaction, user):
            return (
                reaction.message.id == self.message.id
                and str(reaction.emoji) in self.NAV_EMOJIS
                and not user.bot
            )

        while True:
            try:
                reaction, user = await self.bot.wait_for(
                    "reaction_add", timeout=self.timeout, check=check
                )
            except asyncio.TimeoutError:
                try:
                    await self.message.clear_reactions()
                except:
                    pass
                break

            await self.handle_reaction(reaction)
            try:
                await self.message.remove_reaction(reaction.emoji, user)
            except:
                pass

    async def handle_reaction(self, reaction):
        action = self.NAV_EMOJIS[str(reaction.emoji)]
        self.current_page = max(0, min(self.current_page + action, len(self.pages) - 1))
        await self.message.edit(embed=self.pages[self.current_page])
