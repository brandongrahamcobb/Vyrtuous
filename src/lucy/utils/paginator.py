''' paginator.py  The purpose of this program is to override the `menus` Discord pagination implementation.
    Copyright (C) 2024  github.com/brandongrahamcobb

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

class Paginator:
    def __init__(self, bot, ctx, pages):
        self.bot = bot
        self.ctx = ctx
        self.pages = pages
        self.current_page = 0
        self.message = None

    async def start(self):
        if not self.pages:
            await self.ctx.send('There are no tags to display.')
            return
        self.message = await self.ctx.send(embed=self.pages[self.current_page])
        await self.message.add_reaction('⬅️')
        await self.message.add_reaction('➡️')
        await self.message.add_reaction('⏹️')
        def check(reaction, user):
            return (
                user == self.ctx.author
                and reaction.message.id == self.message.id
                and str(reaction.emoji) in ['⬅️', '➡️', '⏹️']
            )
        while True:
            try:
                reaction, user = await self.bot.wait_for('reaction_add', timeout=60.0, check=check)
                if str(reaction.emoji) == '⬅️':
                    if self.current_page > 0:
                        self.current_page -= 1
                        await self.message.edit(embed=self.pages[self.current_page])
                elif str(reaction.emoji) == '➡️':
                    if self.current_page < len(self.pages) - 1:
                        self.current_page += 1
                        await self.message.edit(embed=self.pages[self.current_page])
                elif str(reaction.emoji) == '⏹️':
                    await self.message.clear_reactions()
                    break
                await self.message.remove_reaction(reaction.emoji, user)
            except asyncio.TimeoutError:
                await self.message.clear_reactions()
                break

