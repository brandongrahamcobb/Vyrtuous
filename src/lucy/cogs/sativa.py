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
from lucy.utils.backup import perform_backup, setup_backup_directory
from lucy.utils.helpers import *
from typing import Literal, Optional

import asyncio
import discord

class Sativa(commands.Cog):

    def __init__(self, bot):
        self.bot = bot

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

    @commands.hybrid_command(name='backup', hidden=True)
    @commands.check(is_owner)
    async def backup_task(self, ctx: commands.Context):
        try:
            backup_dir = setup_backup_directory('./backups')
            backup_file = perform_backup(
                db_user='postgres',
                db_name='lucy',
                db_host='localhost',
                backup_dir=backup_dir
            )

            logger.info(f'Backup completed successfully: {backup_file}')
        except Exception as e:
            logger.error(f'Error during database backup: {e}')


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


    @commands.command(name='wipe', hidden=True)
    @commands.check(is_owner)
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
