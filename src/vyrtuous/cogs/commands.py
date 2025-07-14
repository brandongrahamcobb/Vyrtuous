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
from discord import app_commands
import asyncio
import discord
from discord.ext import commands
import inspect
from vyrtuous.utils.handlers.message_service import MessageService, Paginator
from vyrtuous.utils.handlers.predicator import *
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
        self.command_aliases: dict[int, dict[str, dict[str, int]]] = defaultdict(lambda: {'mute': {}, 'unmute': {}})
  
    async def get_available_commands(self, bot, ctx):
        available_commands = []
        for command in bot.commands:
            try:
                if await command.can_run(ctx):
                    available_commands.append(command)
            except commands.CheckFailure:
                continue
        return available_commands
        
    async def cog_load(self):
        await self.load_aliases()

    async def load_aliases(self):
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT guild_id, alias_type, alias_name, channel_id FROM command_aliases')
            for row in rows:
                guild_id = row['guild_id']
                alias_type = row['alias_type']
                alias_name = row['alias_name']
                channel_id = row['channel_id']
                self.command_aliases[guild_id][alias_type][alias_name] = channel_id
                if alias_type == 'mute':
                    cmd = self.create_mute_command(alias_name)
                    self.bot.add_command(cmd)
                elif alias_type == 'unmute':
                    cmd = self.create_unmute_command(alias_name)
                    self.bot.add_command(cmd)

    def create_mute_command(self, command_name: str):
        @commands.hybrid_command(name=command_name, help='Mutes a member in a specific VC.')
        @commands.check(is_moderator)
        async def mute_command(
            ctx,
            member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Optionally include a reason for the mute.')
        ):
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith('<@') and member_input.endswith('>'):
                try:
                    member_id = int(member_input.strip('<@!>'))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
            static_channel_id = self.command_aliases.get(guild_id, {}).get('mute', {}).get(command_name)
            async with self.db_pool.acquire() as conn:
                await conn.execute('''
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
                ''', member_object.id, static_channel_id)
                await conn.execute('''
                    INSERT INTO mute_reasons (guild_id, user_id, reason)
                    VALUES ($1, $2, $3)
                    ON CONFLICT (guild_id, user_id)
                    DO UPDATE SET reason = EXCLUDED.reason
                ''', guild_id, member_object.id, reason)
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await member_object.edit(mute=True)
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been muted in <#{static_channel_id}> with reason {reason}.')
            else:
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been muted in <#{static_channel_id}> with reason {reason}.')
        mute_command.__name__ = f'mute_cmd_{command_name}'
        return mute_command

    @commands.hybrid_command(name='reason', help='Get the reason for a mute or unmute.')
    @commands.check(is_moderator)
    async def reason_command(
        self,
        ctx,
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.')
    ):
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_id, reason
                FROM mute_reasons
                WHERE guild_id = $1 AND user_id = $2
            ''', guild_id, member_object.id)
        if not rows:
            return await self.handler.send_message(ctx, content=f'No mute/unmute history found for {member_object.mention}.')
        lines = []
        for row in rows:
            channel_id = row['channel_id']
            reason = row['reason'] or '*No reason provided*'
            lines.append(f'• <#{channel_id}>: `{reason}`')
        content = f'📄 Mute/Unmute reasons for {member_object.mention}:\n' + '\n'.join(lines)
        await self.handler.send_message(ctx, content=content)

    def create_unmute_command(self, command_name: str):
        @commands.hybrid_command(name=command_name, help='Unmutes a member in a specific VC.')
        @commands.check(is_moderator)
        async def unmute_command(
            ctx,
            member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
            *,
            reason: str = commands.parameter(default='N/A', description='Include a reason for the unmute.')
        ):
            guild_id = ctx.guild.id
            member_id = None
            member_object = None
            if member_input.isdigit():
                member_id = int(member_input)
            elif member_input.startswith('<@') and member_input.endswith('>'):
                try:
                    member_id = int(member_input.strip('<@!>'))
                except ValueError:
                    pass
            if member_id:
                member_object = ctx.guild.get_member(member_id)
            if not member_object:
                return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
            static_channel_id = self.command_aliases.get(guild_id, {}).get('unmute', {}).get(command_name)
            async with self.db_pool.acquire() as conn:
                await conn.execute('''
                    UPDATE users
                    SET mute_channel_ids = array_remove(mute_channel_ids, $2),
                        updated_at = NOW()
                    WHERE user_id = $1
                ''', member_object.id, static_channel_id)
    
                await conn.execute('''
                    INSERT INTO mute_reasons (guild_id, user_id, reason)
                    VALUES ($1, $2, $3)
                    ON CONFLICT (guild_id, user_id)
                    DO UPDATE SET reason = EXCLUDED.reason
                ''', guild_id, member_object.id, reason)
            if member_object.voice and member_object.voice.channel and member_object.voice.channel.id == static_channel_id:
                await member_object.edit(mute=False)
                await self.handler.send_message(ctx, content=f'{member_object.mention} has been unmuted in <#{static_channel_id}>.')
            else:
                await self.handler.send_message(ctx, content=f'{member_object.mention} is no longer marked as muted in <#{static_channel_id}>.')
        unmute_command.__name__ = f'unmute_cmd_{command_name}'
        return unmute_command
    
    # For developers
    @commands.hybrid_command(name='dev', help='Elevates a user's permissions to a bot developer.')
    @commands.check(is_owner)
    async def add_developer(
        self,
        ctx,
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ):
        guild_id = ctx.guild.id
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
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
            ''', member_object.id, guild_id)
        await self.handler.send_message(ctx, content=f'{member_object.mention} has been granted developer rights in this guild.')
        
    @commands.hybrid_command(name='devs', help='Lists all developers.')
    @commands.check(is_owner_or_developer)
    async def list_developers(self, ctx):
        guild = ctx.guild
        pages = []
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT user_id, developer_guild_ids
                FROM users
                WHERE $1 = ANY(developer_guild_ids)
            ''', guild.id)
        if not rows:
            await self.handler.send_message(ctx, content='No developers are configured in this guild.')
            return
        for row in rows:
            user_id = row['user_id']
            user = guild.get_member(user_id)
            name = user.display_name if user else f'User ID {user_id}'
            embed = discord.Embed(
                title=f'Developer: {name}',
                color=discord.Color.blue()
            )
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
        
    @commands.hybrid_command(name='xdev', help='Removes a developer.')
    @commands.check(is_owner_or_developer)
    async def revoke_developer(
        self,
        ctx,
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
    ):
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        guild_id = ctx.guild.id
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member_object.id, guild_id)
        await self.handler.send_message(ctx, content=f'{member_object.mention}'s developer access has been revoked in this guild.')

    # For moderators
    @commands.hybrid_command(name='mod', help='Elevates a user's permission to VC moderator for a specific channel.')
    @commands.check(is_owner_or_developer)
    async def add_moderator(
        self,
        ctx,
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel_input: str = commands.parameter(description='Tag a channel or include its snowflake ID.')
    ):
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        resolved_channel = None
        if channel_input.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel_input))
        elif channel_input.startswith('<#') and channel_input.endswith('>'):
            try:
                channel_id = int(channel_input.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel_input.lower():
                    resolved_channel = vc
                    break
        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='Could not resolve a valid **voice** channel from your input.')
        async with self.db_pool.acquire() as conn:
            await conn.execute('''
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
            ''', member_object.id, resolved_channel.id)
        await self.handler.send_message(ctx, content=f'{member_object.mention} has been granted VC moderator access in {resolved_channel.name}.')

    @commands.hybrid_command(name='mods', help='Lists VC moderators.')
    @commands.check(is_owner_or_developer)
    async def list_moderators(self, ctx):
        guild = ctx.guild
        pages = []
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT user_id, moderator_ids
                FROM users
                WHERE cardinality(moderator_ids) > 0
            ''')
        for row in rows:
            user_id = row['user_id']
            moderator_ids = row['moderator_ids'] or []
            valid_channels = [
                guild.get_channel(cid)
                for cid in moderator_ids
                if (guild.get_channel(cid) and isinstance(guild.get_channel(cid), discord.VoiceChannel))
            ]
            if not valid_channels:
                continue
            user = guild.get_member(user_id)
            display_name = user.display_name if user else f'User ID {user_id}'
            embed = discord.Embed(
                title=f'VC Mod: {display_name}',
                description='\n'.join(f'<#{channel.id}> — {channel.name}' for channel in valid_channels),
                color=discord.Color.blue()
            )
            embed.set_footer(text=f'User ID: {user_id}')
            pages.append(embed)
        if not pages:
            await self.handler.send_message(ctx, content='No VC moderators are configured in this guild.')
            return
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
        
    @commands.hybrid_command(name='xmod', help='Revokes a member's VC moderator role for a given channel.')
    @commands.check(is_owner_or_developer)
    async def revoke_moderator(
        self,
        ctx,
        member_input: str = commands.parameter(description='Tag a user or include their snowflake ID.'),
        channel_input: str = commands.parameter(description='Tag a VC or include its snowflake ID.')
    ):
        member_id = None
        member_object = None
        if member_input.isdigit():
            member_id = int(member_input)
        elif member_input.startswith('<@') and member_input.endswith('>'):
            try:
                member_id = int(member_input.strip('<@!>'))
            except ValueError:
                pass
        if member_id:
            member_object = ctx.guild.get_member(member_id)
        if not member_object:
            return await self.handler.send_message(ctx, content='Could not resolve a valid guild member from your input.')
        resolved_channel = None
        if channel_input.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel_input))
        elif channel_input.startswith('<#') and channel_input.endswith('>'):
            try:
                channel_id = int(channel_input.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel_input.lower():
                    resolved_channel = vc
                    break
        if not isinstance(resolved_channel, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='Could not resolve a valid **voice** channel from your input.')
        async with self.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET moderator_ids = array_remove(moderator_ids, $2),
                    updated_at = NOW()
                WHERE user_id = $1
            ''', member_object.id, resolved_channel.id)
        await self.handler.send_message(ctx, content=f'{member_object.mention} has been revoked moderator access in {resolved_channel.name}.')

    # Aliasing
    @commands.hybrid_command(name='xalias', help='Deletes an alias.')
    @commands.check(is_owner_or_developer)
    async def delete_alias(
        self,
        ctx,
        alias_type: str = commands.parameter(description='Include either `mute` or `unmute`'),
        alias_name: str = commands.parameter(description='Includ an alias name')
    ):
        if alias_type.lower() not in {'mute', 'unmute'}:
            await ctx.send('❌ `alias_type` must be either `mute` or `unmute`.', ephemeral=True)
            return
        if not alias_name.strip():
            await ctx.send('❌ `alias_name` cannot be empty.', ephemeral=True)
            return
        guild_id = ctx.guild.id
        alias_map = self.command_aliases.get(guild_id, {}).get(alias_type.lower(), {})
        if alias_name not in alias_map:
            await ctx.send(f'❌ Alias `{alias_name}` not found in `{alias_type}` for guild `{guild_id}`.', ephemeral=True)
            return
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                'DELETE FROM command_aliases WHERE guild_id = $1 AND alias_type = $2 AND alias_name = $3',
                guild_id, alias_type.lower(), alias_name
            )
        self.command_aliases[guild_id][alias_type.lower()].pop(alias_name, None)
        await self.handler.send_message(ctx, content=f'✅ Deleted alias `{alias_name}` from `{alias_type}`.')
        
    @commands.hybrid_command(name='aliases', help='List all the aliases in the current guild.')
    @commands.check(is_owner_or_developer)
    async def list_aliases(self, ctx):
        guild_id = ctx.guild.id
        aliases = self.command_aliases.get(guild_id, {})
        embed = discord.Embed(title=f'Command Aliases for {ctx.guild.name}')
        for kind in ('mute', 'unmute'):
            lines = [f'`{name}` → <#{cid}>' for name, cid in aliases.get(kind, {}).items()]
            embed.add_field(name=kind.capitalize(), value='\n'.join(lines) or 'None', inline=False)
        await self.handler.send_message(ctx, embed=embed)
        
    @commands.hybrid_command(name='alias', help='Set a mute/unmute alias for a given VC and guild.')
    @commands.check(is_owner_or_developer)
    async def set_alias(
        self,
        ctx,
        alias_type: str = commands.parameter(description='Include either `mute` or `unmute`'),
        alias_name: str = commands.parameter(description='Includ an alias name'),
        channel: str = commands.parameter(description='Include a VC snowflake ID')
    ):
        alias_type = alias_type.lower()
        if alias_type not in {'mute', 'unmute'}:
            await ctx.send('❌ `alias_type` must be either `mute` or `unmute`.', ephemeral=True)
            return
        if not alias_name.strip():
            await ctx.send('❌ `alias_name` cannot be empty.', ephemeral=True)
            return
        resolved_channel = None
        if channel.isdigit():
            resolved_channel = ctx.guild.get_channel(int(channel))
        elif channel.startswith('<#') and channel.endswith('>'):
            try:
                channel_id = int(channel.strip('<#>'))
                resolved_channel = ctx.guild.get_channel(channel_id)
            except ValueError:
                pass
        else:
            for vc in ctx.guild.voice_channels:
                if vc.name.lower() == channel.lower():
                    resolved_channel = vc
                    break
        if not isinstance(resolved_channel, discord.VoiceChannel):
            await self.handler.send_message(ctx, content='❌ Could not resolve a valid VC from your input.')
            return
        guild_id = ctx.guild.id
        async with self.db_pool.acquire() as conn:
            await conn.execute(
                '''
                INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT (guild_id, alias_type, alias_name)
                DO UPDATE SET channel_id = EXCLUDED.channel_id
                ''',
                guild_id, alias_type, alias_name, resolved_channel.id
            )
        self.command_aliases[guild_id][alias_type][alias_name] = resolved_channel.id
        await self.handler.send_message(
            ctx,
            content=f'✅ Alias `{alias_name}` ({alias_type}) set to VC {resolved_channel.mention}.'
        )
        
    @commands.hybrid_command(name='help')
    async def help(
        self,
        ctx,
        *,
        command_name: str = commands.parameter(default=None, description='Include a command name')
    ):
        bot = ctx.bot
        if command_name:
            cmd = bot.get_command(command_name.lower())
            if not cmd:
                await ctx.send(f'❌ Command `{command_name}` not found.')
                return
            if cmd.hidden:
                await ctx.send(f'❌ Command `{command_name}` is hidden.')
                return
            if not await cmd.can_run(ctx):
                await ctx.send(f'❌ You do not have permission to run `{command_name}`.')
                return
            embed = discord.Embed(
                title=f'/{cmd.name}',
                description=cmd.help or 'No description provided.',
                color=discord.Color.blue()
            )
            sig = inspect.signature(cmd.callback)
            parameters = list(sig.parameters.items())[2:]
            for name, param in parameters:
                is_optional = param.kind == inspect.Parameter.KEYWORD_ONLY
                default = param.default
                description = (
                    getattr(default, 'description', None)
                    if isinstance(default, commands.Parameter)
                    else None
                )
                default_value = (
                    getattr(default, 'default', None)
                    if isinstance(default, commands.Parameter)
                    else default
                )
                annotation = (
                    param.annotation.__name__
                    if hasattr(param.annotation, '__name__')
                    else str(param.annotation)
                )
                label = 'Optional' if is_optional else 'Required'
                display = f'Type: `{annotation}`\n{label}'
                if description:
                    display += f'\nDescription: {description}'
                embed.add_field(name=f'`{name}`', value=display, inline=False)
            await ctx.send(embed=embed)
            return
        all_commands = await self.get_available_commands(bot, ctx)
        if not all_commands:
            await ctx.send('❌ No commands available to you.')
            return
        cog_map: dict[str, list[commands.Command]] = {}
        for command in all_commands:
            if command.hidden:
                continue
            cog_map.setdefault(command.cog_name or 'Uncategorized', []).append(command)
        pages = []
        for cog_name in sorted(cog_map):
            commands_in_cog = sorted(cog_map[cog_name], key=lambda c: c.name)
            embed = discord.Embed(title=f'{cog_name} Commands', color=discord.Color.green())
            embed.description = '\n'.join(
                f'**/{cmd.name}** – {cmd.help or 'No description'}' for cmd in commands_in_cog
            )
            pages.append(embed)
        paginator = Paginator(bot, ctx, pages)
        await paginator.start()
        
    
async def setup(bot: commands.Bot):
    cog = Hybrid(bot)
    await bot.add_cog(cog)

