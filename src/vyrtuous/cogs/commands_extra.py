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
from typing import Dict, List, Literal, Optional

import discord
from discord.ext import commands, tasks
from vyrtuous.utils.handlers.message_service import Paginator
from vyrtuous.utils.handlers.predicator import Predicator
from vyrtuous.utils.handlers.sql_manager import perform_backup, setup_backup_directory
from vyrtuous.utils.inc.helpers import *

class Sativa(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.predicator = Predicator(self.bot)

    @staticmethod
    def at_home(bot):
        async def predicate(ctx):
            return ctx.guild is not None and ctx.guild.id == bot.config.get("discord_testing_guild_id")
        return commands.check(predicate)

    @staticmethod
    def is_owner(bot):
        async def predicate(ctx):
            return ctx.guild is not None and (ctx.guild.owner_id == ctx.author.id or ctx.guild.owner_id == bot.config['discord_owner_id'])
        return commands.check(predicate)

    @commands.hybrid_command(name='load', hidden=True)
    @commands.check(at_home)
    async def load(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')


    @commands.hybrid_command(name='reload', hidden=True)
    @commands.check(at_home)
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command(name='sync', hidden=True)
    @commands.check(is_owner)
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


async def setup(bot: commands.bot):
    await bot.add_cog(Sativa(bot))
