''' sativa.py  The purpose of this program is to provide permission-restricted commands to a Discord bot from cd ../../..
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
from discord.ext import commands, tasks
from lucy.utils.helpers import *
from typing import Literal, Optional

import asyncio
import discord

def is_owner():
    async def predicate(ctx):
        return ctx.guild is not None and ctx.guild.owner_id == ctx.author.id
    return commands.check(predicate)

def at_home():
    async def predicate(ctx):
        return ctx.guild is not None and ctx.guild.id == 1300517536001036348
    return commands.check(predicate)

class Sativa(commands.Cog):

    def __init__(self, bot):
        self.bot = bot

    @commands.command(name='sync', hidden=True)
    @is_owner()
    @at_home()
    async def sync(self, ctx: commands.Context, guilds: commands.Greedy[discord.Object], spec: Optional[Literal['~', '*', '^']] = None) -> None:
        if not guilds:
            if spec == '~':
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '*':
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '^':
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
                synced = []
            else:
                synced = await ctx.bot.tree.sync()
            await ctx.send(
                f'Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}'
            )
            return
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        await ctx.send(f'Synced the tree to {ret}/{len(guilds)}.')


    @commands.command(name='wipe', hidden=True)
    @is_owner()
    @at_home()
    async def wipe(self, ctx, option: str = None, limit: int = 100):
        if limit <= 0 or limit > 100:
            return await ctx.send('Limit must be between 1 and 100.')
        
        check_function = None
        if option == 'bot':
            check_function = lambda m: m.author == self.bot.user
        elif option == 'all':
            check_function = lambda m: True  # Allow all messages to be deleted
        elif option == 'user':
            user = ctx.message.mentions[0] if ctx.message.mentions else None
            if user:
                check_function = lambda m: m.author == user
            else:
                return await ctx.send('Please mention a user.')
        elif option == 'commands':
            check_function = lambda m: m.content.startswith(ctx.prefix)
        elif option == 'text':
            await ctx.send('Provide text to delete messages containing it.')
            try:
                msg_text = await self.bot.wait_for('message', timeout=30.0, check=lambda m: m.author == ctx.author)
                check_function = lambda m: msg_text.content in m.content
            except asyncio.TimeoutError:
                return await ctx.send('You took too long to provide text. Cancelling operation.')
        else:
            return await ctx.send('Invalid option.')
        
        total_deleted = 0
        while total_deleted < limit:
            # Purge messages in smaller chunks to avoid hitting rate limits
            deleted = await ctx.channel.purge(limit=min(limit - total_deleted, 10), check=check_function)
            if not deleted:  # Exit loop if no messages were deleted
                break
            total_deleted += len(deleted)
            await asyncio.sleep(1)
        
        if total_deleted > 0:
            await ctx.send(f'Deleted {total_deleted} messages.')
        else:
            await ctx.send('No messages matched the criteria.')

async def setup(bot: commands.bot):
    await bot.add_cog(Sativa(bot))
