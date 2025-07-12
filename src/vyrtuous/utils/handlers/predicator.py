''' predicator.py  The purpose of this program is to provide checks for predication.
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

import discord
from discord.ext import commands
from discord.utils import get
from vyrtuous.utils.inc.setup_logging import logger

class Predicator:
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool

    def at_home(self) -> bool:
        async def predicate(ctx):
            return ctx.guild is not None and ctx.guild.id == self.config['discord_testing_guild_id']
        return commands.check(predicate)

    def release_mode(self) -> bool:
        async def predicate(ctx):
            return (
                ctx.author.id == 154749533429956608 or 
                self.config.get('discord_release_mode', False) or 
                isinstance(ctx.channel, discord.DMChannel)
            )
        return commands.check(predicate)

    async def is_at_home_func(self, guild_id: int) -> bool:
        if guild_id == self.config['discord_testing_guild_id'] or (guild_id in self.config['discord_testing_guild_ids']):
            return True

    def is_developer(self, member: discord.Member) -> bool:
        if member is not None:
            return member.id == self.config['discord_owner_id']
        else:
            return False

    async def is_vegan_user(self, user: discord.User) -> bool:
        guild_ids = self.config['discord_testing_guild_ids']
        for guild_id in guild_ids:
            guild = self.bot.get_guild(guild_id)
            if guild:
                vegan_role = get(guild.roles, name='vegan')
                if vegan_role and vegan_role in user.roles:
                    return True
        return False

    def is_release_mode_func(self, ctx: commands.Context) -> bool:
        return (
            ctx.author.id == 154749533429956608 or
            self.config.get('discord_release_mode', False) or
            isinstance(ctx.channel, discord.DMChannel)
        )
    
    async def is_moderator(self, executor: discord.Member, target: discord.Member) -> bool:
        async with self.db_pool.acquire() as conn:
            row = await conn.fetchrow("""
                SELECT moderator_ids FROM users WHERE user_id = $1
            """, executor.id)

        if not row or not row["moderator_ids"]:
            return False
        return executor.voice.channel.id in row["moderator_ids"]

