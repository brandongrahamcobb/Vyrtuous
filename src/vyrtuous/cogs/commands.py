''' hybrid.py The purpose of this program is to be an extension to a Discord
    bot to provide the command functionality to Py_vyrtuous.
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
from collections import defaultdict
import discord
from discord.ext import commands
from vyrtuous.utils.handlers.message_service import MessageService
from vyrtuous.utils.handlers.predicator import Predicator
from vyrtuous.utils.inc.helpers import *
from types import MethodType

class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.predicator = Predicator(self.bot)
        self.handler = MessageService(self.bot, self.config, self.bot.db_pool)
        self.command_aliases: dict[int, dict[str, dict[str, int]]] = defaultdict(lambda: {"mute": {}, "unmute": {}})
        
    @staticmethod
    def is_owner(bot):
        async def predicate(ctx):
            return ctx.guild is not None and (ctx.guild.owner_id == ctx.author.id or ctx.guild.owner_id == bot.config['discord_owner_id'])
        return commands.check(predicate)
        
    @staticmethod
    def is_owner_or_developer(bot):
        async def predicate(ctx):
            if ctx.guild is None:
                return False
    
            author_id = ctx.author.id
            guild_id = ctx.guild.id
    
            # Guild owner or global owner?
            if ctx.guild.owner_id == author_id or author_id == bot.config["discord_owner_id"]:
                return True
    
            # Check DB if user is a developer in this guild
            async with bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(
                    "SELECT developer_guild_ids FROM users WHERE user_id = $1", author_id
                )
                if row and row["developer_guild_ids"] and guild_id in row["developer_guild_ids"]:
                    return True
    
            return False
    
        return commands.check(predicate)
    
    
    async def load_aliases(self):
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch("SELECT guild_id, alias_type, alias_name, channel_id FROM command_aliases")
            for row in rows:
                guild_id = row["guild_id"]
                alias_type = row["alias_type"]  # 'mute' or 'unmute'
                alias_name = row["alias_name"]
                channel_id = row["channel_id"]
                self.command_aliases[guild_id][alias_type][alias_name] = channel_id

    def create_mute_command(self, command_name: str):
        @commands.command(name=command_name)
        async def mute_command(ctx, member: discord.Member):
            guild_id = ctx.guild.id
            static_channel_id = self.command_aliases.get(guild_id, {}).get("mute", {}).get(command_name)
    
            if not static_channel_id:
                await ctx.send("This mute command is not configured properly.")
                return
    
            if not await self.predicator.is_moderator(ctx.author, static_channel_id):
                await ctx.send("You are not a moderator in this channel.")
                return
    
            async with self.db_pool.acquire() as conn:
                await conn.execute("""
                    INSERT INTO users (user_id, mute_channel_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (user_id) DO UPDATE
                    SET mute_channel_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(u.mute_channel_ids || EXCLUDED.mute_channel_ids)
                        )
                        FROM users u WHERE u.user_id = EXCLUDED.user_id
                    ),
                    updated_at = NOW()
                """, member.id, static_channel_id)
    
            if member.voice and member.voice.channel and member.voice.channel.id == static_channel_id:
                await member.edit(mute=True)
                await ctx.send(f"{member.mention} has been muted in <#{static_channel_id}>.")
            else:
                await ctx.send(f"{member.mention} is now marked to be muted when they join <#{static_channel_id}>.")
    
        mute_command.__name__ = f"mute_cmd_{command_name}"
        return mute_command


    def create_unmute_command(self, command_name: str):
        @commands.command(name=command_name)
        async def unmute_command(ctx, member: discord.Member):
            guild_id = ctx.guild.id
            static_channel_id = self.command_aliases.get(guild_id, {}).get("unmute", {}).get(command_name)
    
            if not static_channel_id:
                await ctx.send("This unmute command is not configured properly.")
                return
    
            if not await self.predicator.is_moderator(ctx.author, static_channel_id):
                await ctx.send("You are not a moderator in this channel.")
                return
    
            async with self.db_pool.acquire() as conn:
                await conn.execute("""
                    UPDATE users
                    SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                """, member.id, static_channel_id)
    
            if member.voice and member.voice.channel and member.voice.channel.id == static_channel_id:
                await member.edit(mute=False)
                await ctx.send(f"{member.mention} has been unmuted in <#{static_channel_id}>.")
            else:
                await ctx.send(f"{member.mention} is no longer marked as muted in <#{static_channel_id}>.")
    
        unmute_command.__name__ = f"unmute_cmd_{command_name}"
        return unmute_command
        
    # For developers
    @commands.command(name="give_dev")
    @commands.check(is_owner)  # Only existing devs/owners can grant
    async def grant_developer(ctx, member: discord.Member):
        guild_id = ctx.guild.id
    
        async with bot.db_pool.acquire() as conn:
            await conn.execute("""
                INSERT INTO users (user_id, developer_guild_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (user_id) DO UPDATE
                SET developer_guild_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(u.developer_guild_ids || EXCLUDED.developer_guild_ids)
                    )
                    FROM users u WHERE u.user_id = EXCLUDED.user_id
                ),
                updated_at = NOW()
            """, member.id, guild_id)
    
        await ctx.send(f"{member.mention} has been granted developer rights in this server.")
        
    @commands.command(name="list_devs")
    @commands.check(is_owner_or_developer)
    async def list_developers(self, ctx):
        guild = ctx.guild
        pages = []
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch("""
                SELECT user_id, developer_guild_ids
                FROM users
                WHERE $1 = ANY(developer_guild_ids)
            """, guild.id)
    
        if not rows:
            await ctx.send("No developers are configured in this server.")
            return
        for row in rows:
            user_id = row["user_id"]
            user = guild.get_member(user_id)
            name = user.display_name if user else f"User ID {user_id}"
    
            embed = discord.Embed(
                title="Developer Access",
                description=f"👤 {name}\n🆔 `{user_id}`",
                color=discord.Color.orange()
            )
            pages.append(embed)
    
        paginator = self.handler.Paginator(self.bot, ctx, pages)
        await paginator.start()
        
    @commands.command(name="revoke_dev")
    @commands.check(is_owner_or_developer)
    async def revoke_developer(ctx, member: discord.Member):
        guild_id = ctx.guild.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute("""
                UPDATE users
                SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            """, member.id, guild_id)
    
        await ctx.send(f"{member.mention}'s developer access has been revoked in this server.")

    # For moderators
    @commands.command(name="give_mod")
    @commands.check(is_owner_or_developer)
    async def add_moderator_channel(self, ctx, member: discord.Member, channel: discord.VoiceChannel):
        async with self.db_pool.acquire() as conn:
            await conn.execute("""
                INSERT INTO users (user_id, moderator_ids)
                VALUES ($1, ARRAY[$2]::BIGINT[])
                ON CONFLICT (user_id) DO UPDATE
                SET moderator_ids = (
                    SELECT ARRAY(
                        SELECT DISTINCT unnest(u.moderator_ids || EXCLUDED.moderator_ids)
                    )
                    FROM users u WHERE u.user_id = EXCLUDED.user_id
                ),
                updated_at = NOW()
            """, member.id, channel.id)
    
        await ctx.send(f"{member.mention} has been granted moderator access in {channel.name}.")

    @commands.command(name="list_mods")
    @commands.check(is_owner_or_developer)
    async def list_moderators(self, ctx):
        guild = ctx.guild
        pages = []
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch("""
                SELECT user_id, moderator_ids
                FROM users
                WHERE cardinality(moderator_ids) > 0
            """)
        for row in rows:
            user_id = row["user_id"]
            moderator_ids = row["moderator_ids"] or []
            valid_channels = [
                guild.get_channel(cid)
                for cid in moderator_ids
                if (guild.get_channel(cid) and isinstance(guild.get_channel(cid), discord.VoiceChannel))
            ]
            if not valid_channels:
                continue
            user = guild.get_member(user_id)
            display_name = user.display_name if user else f"User ID {user_id}"
            embed = discord.Embed(
                title=f"Moderator: {display_name}",
                description="\n".join(f"<#{channel.id}> — {channel.name}" for channel in valid_channels),
                color=discord.Color.blue()
            )
            embed.set_footer(text=f"User ID: {user_id}")
            pages.append(embed)
        if not pages:
            await ctx.send("No moderators are configured in this server.")
            return
        paginator = self.handler.Paginator(self.bot, ctx, pages)
        await paginator.start()

        
    @commands.command(name="revoke_mod")
    @commands.check(is_owner_or_developer)
    async def remove_moderator_channel(self, ctx, member: discord.Member, channel: discord.VoiceChannel):
        async with self.db_pool.acquire() as conn:
            await conn.execute("""
                UPDATE users
                SET moderator_ids = array_remove(moderator_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            """, member.id, channel.id)
        await ctx.send(f"{member.mention} has been revoked moderator access in {channel.name}.")

    # Aliasing
    @commands.command(name="delalias")
    @commands.check(is_owner_or_developer)
    async def delete_alias(self, ctx, alias_type: str, alias_name: str):
        guild_id = ctx.guild.id
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                "DELETE FROM command_aliases WHERE guild_id = $1 AND alias_type = $2 AND alias_name = $3",
                guild_id, alias_type, alias_name
            )
        self.command_aliases[guild_id][alias_type].pop(alias_name, None)
        await ctx.send(f"Deleted alias `{alias_name}` from {alias_type}.")
    
    @commands.command(name="listaliases")
    @commands.check(is_owner_or_developer)
    async def list_aliases(self, ctx):
        guild_id = ctx.guild.id
        aliases = self.command_aliases.get(guild_id, {})
        embed = discord.Embed(title=f"Command Aliases for {ctx.guild.name}")
        for kind in ("mute", "unmute"):
            lines = [f"`{name}` → <#{cid}>" for name, cid in aliases.get(kind, {}).items()]
            embed.add_field(name=kind.capitalize(), value="\n".join(lines) or "None", inline=False)
        await ctx.send(embed=embed)
        
    @commands.command(name="setalias")
    @commands.check(is_owner_or_developer)
    async def set_alias(self, ctx, alias_type: str, alias_name: str, channel: discord.TextChannel):
        """Set or update an alias for mute/unmute."""
        if alias_type not in ("mute", "unmute"):
            return await ctx.send("Alias type must be `mute` or `unmute`.")
    
        guild_id = ctx.guild.id
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                """
                INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_id, alias_type, alias_name)
                DO UPDATE SET channel_id = EXCLUDED.channel_id
                """,
                guild_id, alias_type, alias_name, channel.id
            )
    
        self.command_aliases[guild_id][alias_type][alias_name] = channel.id
        await ctx.send(f"Alias `{alias_name}` ({alias_type}) set to channel {channel.mention}.")


    @commands.command(name="help")
    @commands.check(is_owner_or_developer)  # apply the check here
    async def restricted_help(self, ctx):
        """Displays help menu, but only for owners or developers."""
        embed = discord.Embed(title="Py_vyrtuous Command Help", color=discord.Color.blurple())
    
        embed.add_field(name="give_mod", value="Grants a user mod rights in a voice channel.", inline=False)
        embed.add_field(name="revoke_mod", value="Revokes a user's mod rights in a voice channel.", inline=False)
        embed.add_field(name="give_dev", value="Grants a user developer rights in this server.", inline=False)
        embed.add_field(name="remove_dev", value="Revokes a user's developer rights.", inline=False)
        embed.add_field(name="setalias", value="Set an alias for a mute/unmute command.", inline=False)
        embed.add_field(name="delalias", value="Delete an alias for a mute/unmute command.", inline=False)
        embed.add_field(name="listaliases", value="List all aliases configured in this server.", inline=False)
    
        await ctx.send(embed=embed)
        

async def setup(bot: commands.Bot):
    cog = Hybrid(bot)
    await cog.load_aliases()
    await bot.add_cog(Hybrid(bot))
    
    @bot.check
    async def restrict_help(ctx):
        if ctx.command.name == "help":
            author_id = ctx.author.id
            guild_id = ctx.guild.id if ctx.guild else None
    
            if guild_id and (
                ctx.guild.owner_id == author_id or author_id == bot.config["discord_owner_id"]
            ):
                return True
    
            async with bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(
                    "SELECT developer_guild_ids FROM users WHERE user_id = $1", author_id
                )
                return row and guild_id in (row["developer_guild_ids"] or [])
    
        return True
