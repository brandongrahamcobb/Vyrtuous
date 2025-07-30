''' message.py  The purpose of this program is to handle messages in Discord.

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
import os

import discord
from discord.ext import commands
from vyrtuous.inc.helpers import *

os.makedirs(DIR_TEMP, exist_ok=True)

class DiscordMessageService:
    def __init__(self, bot, db_pool):
        self.bot = bot

    async def send_dm(self, ctx: commands.Context, *, content: str = None, file: discord.File = None, embed: discord.Embed = None):
        member = ctx.author
        channel = await member.create_dm()
        await self._send_message(channel.send, content=content, file=file, embed=embed)

    async def send_message(self, ctx: commands.Context, *, content: str = None, file: discord.File = None, embed: discord.Embed = None):
        can_send = (
            ctx.guild
            and isinstance(ctx.channel, discord.abc.GuildChannel)
            and ctx.channel.permissions_for(ctx.guild.me).send_messages
        )
        if can_send:
            try:
                await self._send_message(lambda **kwargs: ctx.reply(**kwargs), content=content, file=file, embed=embed)
            except discord.HTTPException as e:
                if e.code == 50035:  # Invalid Form Body due to message_reference
                    await self._send_message(lambda **kwargs: ctx.send(**kwargs), content=content, file=file, embed=embed)
                else:
                    raise
        else:
            await self.send_dm(ctx.author, content=content, file=file, embed=embed)


    async def _send_message(self, send_func, *, content: str = None, file: discord.File = None, embed: discord.Embed = None):
        kwargs = {}
        if content:
            kwargs['content'] = content
        if file:
            kwargs['file'] = file
        if embed:
            kwargs['embed'] = embed
        await send_func(**kwargs)

class Paginator:
    def __init__(self, bot, ctx, pages):
        self.bot = bot
        self.ctx = ctx
        self.pages = pages
        self.current_page = 0
        self.message = None

    async def start(self):
        if not self.pages:
            await self.ctx.send('There are no pages to display.')
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
