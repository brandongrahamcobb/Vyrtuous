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
import logging

logger = logging.getLogger(__name__)
class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.db_pool = bot.db_pool
        self.predicator = Predicator(self.bot)
        self.handler = MessageService(self.bot, self.config, self.bot.db_pool)
        self.command_aliases: dict[int, dict[str, dict[str, int]]] = defaultdict(lambda: {"mute": {}, "unmute": {}})

    async def cog_load(self):
        await self.load_aliases()
    
    @staticmethod
    def is_moderator(bot):
        async def predicate(ctx):
            if not ctx.guild or not ctx.author or not isinstance(ctx.channel, discord.VoiceChannel):
                return False
            channel_id = ctx.channel.id
            async with ctx.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow("""
                    SELECT moderator_ids FROM users WHERE user_id = $1
                """, ctx.author.id)
            if not row or not row["moderator_ids"]:
                return False
            return channel_id in row["moderator_ids"]
        return commands.check(predicate)
    
    
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
            if ctx.guild.owner_id == author_id or author_id == bot.config["discord_owner_id"]:
                return True
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
                alias_type = row["alias_type"]
                alias_name = row["alias_name"]
                channel_id = row["channel_id"]
                self.command_aliases[guild_id][alias_type][alias_name] = channel_id
                if alias_type == "mute":
                    cmd = self.create_mute_command(alias_name)
                    self.bot.add_command(cmd)
                elif alias_type == "unmute":
                    cmd = self.create_unmute_command(alias_name)
                    self.bot.add_command(cmd)

    def create_mute_command(self, command_name: str):
        @commands.command(name=command_name)
        async def mute_command(ctx, member_input: str, *, reason: str = None):
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith("<@") and member_input.endswith(">"):
                try:
                    member_id = int(member_input.strip("<@!>"))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content="Could not resolve a valid guild member from your input.")
            static_channel_id = self.command_aliases.get(guild_id, {}).get("mute", {}).get(command_name)
            if not static_channel_id:
                await self.handler.send_message(ctx, content="This mute command is not configured properly.")
                return
            if not await self.predicator.is_moderator(ctx.author, static_channel_id):
                await self.handler.send_message(ctx, content="You are not a moderator in this channel.")
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
                """, member_object.id, static_channel_id)
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await member_object.edit(mute=False)
                await self.handler.send_message(ctx, content=f"{member_object.mention} has been muted in <#{static_channel_id}>.")
            else:
                await self.handler.send_message(ctx, content=f"{member_object.mention} is marked as muted in <#{static_channel_id}>.")
        mute_command.__name__ = f"mute_cmd_{command_name}"
        return mute_command


    def create_reason_command(self, command_name: str):
        @commands.command(name=command_name)
        async def reason_command(ctx, member_input: str):
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
    
            # Resolve user input
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith("<@") and member_input.endswith(">"):
                try:
                    member_id = int(member_input.strip("<@!>"))
                except ValueError:
                    pass
    
            # Resolve member object
            if member_id:
                member_object = ctx.guild.get_member(member_id)
    
            if not member_object:
                return await self.handler.send_message(ctx, content="Could not resolve a valid guild member from your input.")
    
            # Get static_channel_id associated with the alias command (e.g., "mute" or     "unmute")
            static_channel_id = self.command_aliases.get(guild_id, {}).get("mute", {}).get(command_name)
            if not static_channel_id:
                return await self.handler.send_message(ctx, content="This mute command is not configured properly.")
    
            # Check moderator permissions
            if not await self.predicator.is_moderator(ctx.author, static_channel_id):
                return await self.handler.send_message(ctx, content="You are not a moderator in this channel.")
    
            # Fetch the mute reason from the database
            async with self.db_pool.acquire() as conn:
                row = await conn.fetchrow("""
                    SELECT reason FROM mute_reasons
                    WHERE guild_id = $1 AND user_id = $2
                """, guild_id, member_object.id)
    
            if row and row["reason"]:
                await self.handler.send_message(ctx, content=f"📝 Mute reason for {member_object.mention}: `{row['reason']}`")
            else:
                await self.handler.send_message(ctx, content=f"No mute reason found for {member_object.mention}.")
    
        return reason_command

    def create_unmute_command(self, command_name: str):
        @commands.command(name=command_name)
        async def unmute_command(ctx, member_input: str, *, reason: str = None):
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith("<@") and member_input.endswith(">"):
                try:
                    member_id = int(member_input.strip("<@!>"))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content="Could not resolve a valid guild member from your input.")
            static_channel_id = self.command_aliases.get(guild_id, {}).get("unmute", {}).get(command_name)
            if not static_channel_id:
                await self.handler.send_message(ctx, content="This unmute command is not configured properly.")
                return
            if not await self.predicator.is_moderator(ctx.author, static_channel_id):
                await self.handler.send_message(ctx, content="You are not a moderator in this channel.")
                return
            async with self.db_pool.acquire() as conn:
                await conn.execute("""
                    UPDATE users
                    SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                """, member_object.id, static_channel_id)
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await member_object.edit(mute=False)
                await self.handler.send_message(ctx, content=f"{member_object.mention} has been unmuted in <#{static_channel_id}>.")
            else:
                await self.handler.send_message(ctx, content=f"{member_object.mention} is no longer marked as muted in <#{static_channel_id}>.")
    
        unmute_command.__name__ = f"unmute_cmd_{command_name}"
        return unmute_command
        
    # For developers
    @commands.command(name="give_dev")
    @commands.check(is_owner)  # Only existing devs/owners can grant
    async def grant_developer(ctx, member_input: str):
        guild_id = ctx.guild.id
        member_id = None
        member_object = None

        # Try resolving by raw ID
        if member_input.isdigit():
            member_id = int(member_input)
        # Try resolving by mention
        elif member_input.startswith("<@") and member_input.endswith(">"):
            try:
                member_id = int(member_input.strip("<@!>"))
            except ValueError:
                pass

        # If member_id was successfully parsed, try to get the member object
        if member_id:
            member_object = ctx.guild.get_member(member_id)

        if not member_object:
            return await self.handler.send_message(ctx, content="Could not resolve a valid guild member from your input.")

    
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
            """, member_object.id, guild_id)
    
        await self.handler.send_message(ctx, content=f"{member_object.mention} has been granted developer rights in this server.")
        
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
            await self.handler.send_message(ctx, content="No developers are configured in this server.")
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
    async def revoke_developer(ctx, member_input: str):
        member_id = None
        member_object = None

        # Try resolving by raw ID
        if member_input.isdigit():
            member_id = int(member_input)
        # Try resolving by mention
        elif member_input.startswith("<@") and member_input.endswith(">"):
            try:
                member_id = int(member_input.strip("<@!>"))
            except ValueError:
                pass

        # If member_id was successfully parsed, try to get the member object
        if member_id:
            member_object = ctx.guild.get_member(member_id)

        if not member_object:
            return await self.handler.send_message(ctx, content="Could not resolve a valid guild member from your input.")

        guild_id = ctx.guild.id
        async with bot.db_pool.acquire() as conn:
            await conn.execute("""
                UPDATE users
                SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            """, member_object.id, guild_id)
    
        await self.handler.send_message(ctx, content=f"{member_object.mention}'s developer access has been revoked in this server.")

    # For moderators
    @commands.command(name="give_mod")
    @commands.check(is_owner_or_developer)
    async def add_moderator_channel(self, ctx, member_input: str, channel_input: str):
        """Add a voice channel to the user's moderator_ids array."""
        member_id = None
        member_object = None

        # Try resolving by raw ID
        if member_input.isdigit():
            member_id = int(member_input)
        # Try resolving by mention
        elif member_input.startswith("<@") and member_input.endswith(">"):
            try:
                member_id = int(member_input.strip("<@!>"))
            except ValueError:
                pass

        # If member_id was successfully parsed, try to get the member object
        if member_id:
            member_object = ctx.guild.get_member(member_id)

        if not member_object:
            return await self.handler.send_message(ctx, content="Could not resolve a valid guild member from your input.")

        # Resolve string to VoiceChannel
        resolved_channel = None
    
        if channel_input.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel_input))
        elif channel_input.startswith("<#") and channel_input.endswith(">"):
            try:
                channel_id = int(channel_input.strip("<#>"))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel_input.lower():
                    resolved_channel = vc
                    break
    
        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content="Could not resolve a valid **voice** channel from yourinput.")

    # DB insert/update
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
            """, member_object.id, resolved_channel.id)

        await self.handler.send_message(ctx, content=f"{member_object.mention} has been granted moderator access in {resolved_channel.name}.")

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
            await self.handler.send_message(ctx, content="No moderators are configured in this server.")
            return
        paginator = self.handler.Paginator(self.bot, ctx, pages)
        await paginator.start()

        
    @commands.command(name="revoke_mod")
    @commands.check(is_owner_or_developer)
    async def remove_moderator_channel(self, ctx, member_input: str, channel_input: str):
        """Remove a voice channel from the user's moderator list."""
        member_id = None
        member_object = None

        # Try resolving by raw ID
        if member_input.isdigit():
            member_id = int(member_input)
        # Try resolving by mention
        elif member_input.startswith("<@") and member_input.endswith(">"):
            try:
                member_id = int(member_input.strip("<@!>"))
            except ValueError:
                pass

        # If member_id was successfully parsed, try to get the member object
        if member_id:
            member_object = ctx.guild.get_member(member_id)

        if not member_object:
            return await self.handler.send_message(ctx, content="Could not resolve a valid guild member from your input.")

        # Resolve channel from string input
        resolved_channel = None
    
        if channel_input.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel_input))
        elif channel_input.startswith("<#") and channel_input.endswith(">"):
            try:
                channel_id = int(channel_input.strip("<#>"))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel_input.lower():
                    resolved_channel = vc
                    break
        

        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content="Could not resolve a valid **voice** channel from your input.")

        async with self.db_pool.acquire() as conn:
            await conn.execute("""
                UPDATE users
                SET moderator_ids = array_remove(moderator_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            """, member_object.id, resolved_channel.id)

        await self.handler.send_message(ctx, content=f"{member_object.mention} has been revoked moderator access in {resolved_channel.name}.")

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
        await self.handler.send_message(ctx, content=f"Deleted alias `{alias_name}` from {alias_type}.")
    
    @commands.command(name="listaliases")
    @commands.check(is_owner_or_developer)
    async def list_aliases(self, ctx):
        guild_id = ctx.guild.id
        aliases = self.command_aliases.get(guild_id, {})
        embed = discord.Embed(title=f"Command Aliases for {ctx.guild.name}")
        for kind in ("mute", "unmute"):
            lines = [f"`{name}` → <#{cid}>" for name, cid in aliases.get(kind, {}).items()]
            embed.add_field(name=kind.capitalize(), value="\n".join(lines) or "None", inline=False)
        await self.handler.send_message(ctx, embed=embed)
        
    @commands.command(name="setalias")
    @commands.check(is_owner_or_developer)
    async def set_alias(self, ctx, alias_type: str, alias_name: str, channel: str):
        """Set or update an alias for mute/unmute."""
    
        if alias_type not in ("mute", "unmute"):
            return await self.handler.send_message(ctx, content="Alias type must be `mute` or `unmute`.")
    
        # Convert str to VoiceChannel object
        resolved_channel = None
    
        if channel.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel))
        elif channel.startswith("<#") and channel.endswith(">"):
            try:
                channel_id = int(channel.strip("<#>"))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel.lower():
                    resolved_channel = vc
                    break
    
        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content="Could not resolve a valid **voice** channel from your input.")

        guild_id = ctx.guild.id
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                """
                INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_id, alias_type, alias_name)
                DO UPDATE SET channel_id = EXCLUDED.channel_id
                """,
                guild_id, alias_type, alias_name, resolved_channel.id
            )

        self.command_aliases[guild_id][alias_type][alias_name] = resolved_channel.id
        await self.handler.send_message(ctx, content=f"Alias `{alias_name}` ({alias_type}) set to voice channel {resolved_channel.mention}.")

async def setup(bot: commands.Bot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)
    
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
