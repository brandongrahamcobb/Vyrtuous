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
from typing import Literal, Optional

import utils.helpers as helpers
import discord

def is_owner():
    async def predicate(ctx):
        return ctx.guild is not None and ctx.guild.owner_id == ctx.author.id
    return commands.check(predicate)

class Sativa(commands.Cog):

    def __init__(self, bot):
        self.bot = bot

    @commands.command(name='sync', hidden=True)
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
    async def wipe(self, ctx, option=None, limit: int = 100):
        if limit <= 0 or limit > 100:
            await ctx.send('Please provide a limit between 1 and 100.')
            return
        if option == 'bot':
            def is_bot_message(message):
                return message.author == self.bot.user
            deleted = await self.purge_messages(ctx, limit, is_bot_message)
            await ctx.send(f'Deleted {deleted} bot messages.')
        elif option == 'all':
            deleted = await ctx.channel.purge(limit=limit)
            await ctx.send(f'Deleted all messages.')
        else:
            await ctx.send('Invalid option. Use `bot`, `all`, `user`, or `commands`.')

async def setup(bot: commands.bot):
    await bot.add_cog(Sativa(bot))
