''' check_service.py

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

from discord.ext import commands
from discord import app_commands
from typing import Optional, Tuple
from vyrtuous.config import Config
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.setup_logging import logger
import discord
config = Config.get_config()
      
class NotCoordinator(commands.CheckFailure):
    def __init__(self, message='You are not a coordinator in the requested channel.'):
        super().__init__(message)

class NoInheritedPermissionsForRole(commands.CheckFailure):
    def __init__(self, message='You are not a coordinator in the requested channel.'):
        super().__init__(message)

class NotModerator(commands.CheckFailure):
    def __init__(self, message='You are not a moderator in the requested channel.'):
        super().__init__(message)

class NotGuildOwner(commands.CheckFailure):
    def __init__(self, message='You are not the guild owner and cannot do this.'):
        super().__init__(message)

class CantChoseGuildOwner(commands.CheckFailure):
    def __init__(self, message='You selected the guild owner and cannot do this.'):
        super().__init__(message)

class NotInAGuild(commands.CheckFailure):
    def __init__(self, message='You did not send your message in a guild and cannot do this.'):
        super().__init__(message)

class NotAValidDatabaseQuery(commands.CheckFailure):
    def __init__(self, message='Database query is invalid.'):
        super().__init__(message)

class NoPermissionsForRoleSetup(commands.CheckFailure):
    def __init__(self, message='No permissions for this role are setup.'):
        super().__init__(message)

class NotSystemOwner(commands.CheckFailure):
    def __init__(self, message='You are not the bot owner.'):
        super().__init__(message)

class NotDeveloper(commands.CheckFailure):
    def __init__(self, message='You are not a developer in this guild.'):
        super().__init__(message)

class NoCommandAlias(commands.CheckFailure):
    def __init__(self, message='This command alias is not mapped to a target channel.'):
        super().__init__(message)
                                    
class NotAtHome(commands.CheckFailure):
    def __init__(self, message='You are not in the home Discord!.'):
        super().__init__(message)

async def at_home(ctx) -> bool:
    if ctx.guild is not None and ctx.guild.id == int(config['discord_testing_guild_id']):
        return True
    raise NotAtHome()

async def has_command_alias(ctx) -> bool:
    if not ctx.guild or not ctx.command:
        raise NoCommandAlias('Command must be used in a guild and must be a valid command.')
    bot = ctx.bot
    guild_id = ctx.guild.id
    command_name = ctx.command.name.lower()
    for alias_type in ('role','unrole','mute','unmute','ban','unban','flag','unflag','tmute','untmute','cow','uncow'):
        alias_map = bot.command_aliases.get(guild_id, {}).get(alias_type, {})
        if command_name in alias_map:
            ctx._target_channel_id = alias_map[command_name]
            ctx._target_room_name = None
            return True
        temp_map = bot.command_aliases.get(guild_id, {}).get('temp_room_aliases', {}).get(alias_type, {})
        if command_name in temp_map:
            ctx._target_room_name = temp_map[command_name].get('room_name')
            ctx._target_channel_id = None
            return True
    raise NoCommandAlias(f'Command alias `{command_name}` is not mapped to any target channel.')


async def is_guild_owner_block(ctx, user_id: int):
    if user_id == ctx.guild.owner.id:
        raise CantChoseGuildOwner()
    return True

async def is_owner_block(ctx, user_id: int):
    errors = []
    for check in (is_guild_owner_block):
        try:
            await check(ctx, user_id)
            return True
        except commands.CheckFailure as e:
            errors.append(str(e))
    raise commands.CheckFailure('\n'.join(f'❌ {msg}' for msg in errors))
    
async def send_check_failure_embed(ctx: commands.Context, error: commands.CheckFailure, *, title: str = 'Permission Check Failed'):
    embed = discord.Embed(
        title=title,
        description=error.args[0] if error.args else 'An unknown error occurred.',
        color=discord.Color.red()
    )
    embed.set_footer(text='Contact a bot admin if you think this is a mistake.')
    await ctx.send(embed=embed)

async def is_developer(ctx):
    if ctx.guild is None:
        raise NotDeveloper('Command must be used in a guild.')
    async with ctx.bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT developer_guild_ids FROM users WHERE discord_snowflake = $1', ctx.author.id
        )
    if not row or not row.get('developer_guild_ids') or ctx.guild.id not in row.get('developer_guild_ids', []):
        raise NotDeveloper('You are not a developer in this guild.')
    return True
    
async def is_developer_member(member: discord.Member, bot: commands.Bot) -> bool:
    async with bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT developer_guild_ids FROM users WHERE discord_snowflake = $1', member.id
        )
    if not row or not row.get('developer_guild_ids') or member.guild.id not in row.get('developer_guild_ids', []):
        return False
    return True

async def is_guild_owner(ctx):
    if ctx.guild is None:
        raise NotInAGuild('Command must be used in a guild.')
    if ctx.guild.owner_id != ctx.author.id:
        raise NotGuildOwner()
    return True

async def is_system_owner(ctx):
    system_owner_id = int(ctx.bot.config['discord_owner_id'])
    if ctx.author.id != system_owner_id:
        raise NotSystemOwner()
    return True
    
async def is_owner(ctx: commands.Context, member_id: int):
    system_owner_id = int(ctx.bot.config['discord_owner_id'])
    if ctx.guild.owner_id == member_id or system_owner_id == member_id:
        return True
    return False

async def is_guild_owner_member(member: discord.Member) -> bool:
    if member.guild is None:
        raise NotInAGuild('Member is not in a guild.')
    if member.guild.owner_id != member.id:
        return False
    return True

async def is_system_owner_member(member: discord.Member, bot: commands.Bot) -> bool:
    system_owner_id = int(bot.config['discord_owner_id'])
    if member.id != system_owner_id:
        return False
    return True

async def is_owner_member(member: discord.Member, bot: commands.Bot) -> bool:
    system_owner_id = int(bot.config['discord_owner_id'])
    if member.guild is not None and member.guild.owner_id == member.id:
        return True
    if member.id == system_owner_id:
        return True
    return False

def is_owner_predicator():
    async def predicate(ctx: commands.Context):
        for check in (is_system_owner, is_guild_owner):
            try:
                if await check(ctx):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not system owner or guild owner.')
    predicate._permission_level = 'Owner'
    return commands.check(predicate)

def is_owner_developer_predicator():
    async def predicate(ctx: commands.Context):
        for check in (is_system_owner, is_guild_owner, is_developer):
            try:
                if await check(ctx):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not system owner, guild owner, or developer.')
    predicate._permission_level = 'Developer'
    return commands.check(predicate)

def is_owner_developer_coordinator_moderator_predicator(alias_type: Optional[str] = None):
    async def predicate(ctx: commands.Context):
        guild_aliases = getattr(ctx.bot, 'command_aliases', {}).get(ctx.guild.id, {}) if ctx.guild else {}
        target_channel_id = None
        target_room_name = None
        if alias_type and ctx.guild:
            alias_data_role = guild_aliases.get('role_aliases', {}).get(alias_type, {})
            alias_data_chan = guild_aliases.get('channel_aliases', {}).get(alias_type, {})
            for alias_data in (alias_data_role, alias_data_chan):
                if alias_data:
                    target_channel_id = alias_data.get('channel_id')
                    break
            temp_aliases = guild_aliases.get('temp_room_aliases', {}).get(alias_type, {})
            if temp_aliases:
                alias_name = ctx.command.name.lower()
                if alias_name in temp_aliases:
                    target_room_name = temp_aliases[alias_name].get('room_name')
            if not target_channel_id and not target_room_name and ctx.command:
                async with ctx.bot.db_pool.acquire() as conn:
                    row = await conn.fetchrow('SELECT channel_id FROM command_aliases WHERE guild_id=$1 AND alias_name=$2', ctx.guild.id, ctx.command.name.lower())
                if row:
                    target_channel_id = row['channel_id']
        async with ctx.bot.db_pool.acquire() as conn:
            user_row = await conn.fetchrow('SELECT coordinator_channel_ids,moderator_channel_ids,coordinator_room_names,moderator_room_names FROM users WHERE discord_snowflake=$1', ctx.author.id)
        if user_row:
            c_chan = user_row.get('coordinator_channel_ids') or []
            m_chan = user_row.get('moderator_channel_ids') or []
            c_room = user_row.get('coordinator_room_names') or []
            m_room = user_row.get('moderator_room_names') or []
            if alias_type is None:
                if c_chan or m_chan or c_room or m_room:
                    return True
            else:
                if target_channel_id in c_chan or target_channel_id in m_chan:
                    return True
                if target_room_name in c_room or target_room_name in m_room:
                    return True
        for check in (is_system_owner, is_guild_owner, is_developer):
            try:
                if await check(ctx): return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not system owner, guild owner, developer, coordinator, or moderator.')
    predicate._permission_level = 'Moderator'
    return commands.check(predicate)
    
def is_owner_developer_coordinator_moderator_role_predicator(alias_type: Optional[str] = None):
    async def predicate(ctx: commands.Context):
        db_pool=ctx.bot.db_pool
        guild_aliases=getattr(ctx.bot,'command_aliases',{}).get(ctx.guild.id,{}) if ctx.guild else {}
        target_channel_id=None
        if alias_type and ctx.guild:
            alias_data_chan=guild_aliases.get('channel_aliases',{}).get(alias_type,{})
            alias_data_role=guild_aliases.get('role_aliases',{}).get(alias_type,{})
            for alias_data in (alias_data_role,alias_data_chan):
                if alias_data:
                    target_channel_id=alias_data.get('channel_id')
                    break
            if not target_channel_id and ctx.command:
                async with db_pool.acquire() as conn:
                    row=await conn.fetchrow('SELECT channel_id FROM command_aliases WHERE guild_id=$1 AND alias_name=$2',ctx.guild.id,ctx.command.name.lower())
                if row: target_channel_id=row['channel_id']
        async with db_pool.acquire() as conn:
            role_rows=await conn.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
            user_row=await conn.fetchrow('SELECT coordinator_channel_ids,moderator_channel_ids FROM users WHERE discord_snowflake=$1',ctx.author.id)
        role_ids=[r['role_id'] for r in role_rows] if role_rows else []
        user_role_ids=[role.id for role in getattr(ctx.author,'roles',[])]
        if user_row:
            if alias_type is None:
                if user_row.get('coordinator_channel_ids') or user_row.get('moderator_channel_ids'): return True
            else:
                if target_channel_id in (user_row.get('coordinator_channel_ids') or []): return True
                if target_channel_id in (user_row.get('moderator_channel_ids') or []): return True
        if any(rid in role_ids for rid in user_role_ids): return True
        for check in (is_system_owner,is_guild_owner,is_developer):
            try:
                if await check(ctx): return True
            except commands.CheckFailure: continue
        raise commands.CheckFailure('You are not system owner, guild owner, developer, coordinator, moderator, or team member.')
    predicate._permission_level='Moderator'
    return commands.check(predicate)

def is_owner_developer_coordinator_predicator(alias_type: Optional[str] = None):
    async def predicate(ctx: commands.Context):
        guild_aliases = getattr(ctx.bot, 'command_aliases', {}).get(ctx.guild.id, {}) if ctx.guild else {}
        target_channel_id = None
        target_room_name = None
        if alias_type and ctx.guild:
            alias_data_role = guild_aliases.get('role_aliases', {}).get(alias_type, {})
            alias_data_chan = guild_aliases.get('channel_aliases', {}).get(alias_type, {})
            for alias_data in (alias_data_role, alias_data_chan):
                if alias_data:
                    target_channel_id = alias_data.get('channel_id')
                    break
            temp_aliases = guild_aliases.get('temp_room_aliases', {}).get(alias_type, {})
            if temp_aliases:
                alias_name = ctx.command.name.lower()
                if alias_name in temp_aliases:
                    target_room_name = temp_aliases[alias_name].get('room_name')
            if not target_channel_id and not target_room_name and ctx.command:
                async with ctx.bot.db_pool.acquire() as conn:
                    row = await conn.fetchrow(
                        'SELECT channel_id FROM command_aliases WHERE guild_id=$1 AND alias_name=$2',
                        ctx.guild.id, ctx.command.name.lower()
                    )
                if row:
                    target_channel_id = row['channel_id']
        async with ctx.bot.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake=$1',
                ctx.author.id
            )
        if user_row:
            c_chan = user_row.get('coordinator_channel_ids') or []
            c_room = user_row.get('coordinator_room_names') or []
            if alias_type is None:
                if c_chan or c_room:
                    return True
            else:
                if target_channel_id in c_chan:
                    return True
                if target_room_name in c_room:
                    return True
        for check in (is_system_owner, is_guild_owner, is_developer):
            try:
                if await check(ctx): return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not system owner, guild owner, developer, or coordinator.')
    predicate._permission_level = 'Coordinator'
    return commands.check(predicate)

async def check_owner_dev_coord_mod(ctx: commands.Context, channel: Optional[discord.abc.GuildChannel] = None) -> Tuple[bool, bool]:
    is_owner_or_dev = False
    is_mod_or_coord = False
    guild = ctx.guild
    guild_id = guild.id if guild else None
    user_id = ctx.author.id
    target_channel_id = getattr(channel, 'id', None) or getattr(ctx, '_target_channel_id', None) or getattr(ctx.channel, 'id', None)
    target_room_name = getattr(channel, 'name', None)
    try:
        if await is_system_owner(ctx): is_owner_or_dev = True
    except commands.CheckFailure:
        try:
            if guild and await is_guild_owner(ctx): is_owner_or_dev = True
        except commands.CheckFailure:
            try:
                if guild and await is_developer(ctx): is_owner_or_dev = True
            except commands.CheckFailure:
                pass
    if guild:
        async with ctx.bot.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids,moderator_channel_ids,coordinator_room_names,moderator_room_names FROM users WHERE discord_snowflake=$1',
                user_id
            )
        if user_row:
            c_chan = user_row.get('coordinator_channel_ids') or []
            m_chan = user_row.get('moderator_channel_ids') or []
            c_room = user_row.get('coordinator_room_names') or []
            m_room = user_row.get('moderator_room_names') or []
            if target_channel_id in c_chan or target_channel_id in m_chan:
                is_mod_or_coord = True
            elif target_room_name in c_room or target_room_name in m_room:
                is_mod_or_coord = True
    return is_owner_or_dev, is_mod_or_coord
    
async def check_owner_dev_coord(ctx: commands.Context, channel: Optional[discord.abc.GuildChannel] = None) -> Tuple[bool, bool]:
    is_owner_or_dev = False
    is_coord = False
    guild = ctx.guild
    user_id = ctx.author.id
    target_channel_id = getattr(channel, 'id', None) or getattr(ctx, '_target_channel_id', None) or getattr(ctx.channel, 'id', None)
    target_room_name = getattr(channel, 'name', None)
    try:
        if await is_system_owner(ctx): is_owner_or_dev = True
    except commands.CheckFailure:
        try:
            if guild and await is_guild_owner(ctx): is_owner_or_dev = True
        except commands.CheckFailure:
            try:
                if guild and await is_developer(ctx): is_owner_or_dev = True
            except commands.CheckFailure:
                pass
    if guild:
        async with ctx.bot.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids,coordinator_room_names FROM users WHERE discord_snowflake=$1',
                user_id
            )
        if user_row:
            c_chan = user_row.get('coordinator_channel_ids') or []
            c_room = user_row.get('coordinator_room_names') or []
            if target_channel_id in c_chan or target_room_name in c_room:
                is_coord = True
    return is_owner_or_dev, is_coord
    
async def check_owner_dev_coord_app(
    interaction: discord.Interaction,
    channel: Optional[discord.abc.GuildChannel] = None
) -> Tuple[bool, bool]:
    is_owner_or_dev = False
    is_coord = False
    guild = interaction.guild
    user_id = interaction.user.id
    target_channel_id = getattr(channel, 'id', None) or getattr(interaction, '_target_channel_id', None) or getattr(interaction.channel, 'id', None)
    target_room_name = getattr(channel, 'name', None)
    try:
        if await is_system_owner_app(interaction):
            is_owner_or_dev = True
    except Exception:
        try:
            if guild and await is_guild_owner_app(interaction):
                is_owner_or_dev = True
        except Exception:
            try:
                if guild and await is_developer_app(interaction):
                    is_owner_or_dev = True
            except Exception:
                pass
    if guild:
        async with interaction.client.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake=$1',
                user_id
            )
        if user_row:
            c_chan = user_row.get('coordinator_channel_ids') or []
            c_room = user_row.get('coordinator_room_names') or []
            if target_channel_id in c_chan or target_room_name in c_room:
                is_coord = True
    return is_owner_or_dev, is_coord

async def is_coordinator(ctx, target_channel_id: int):
    if ctx.guild is None:
        raise NotCoordinator('Command must be used in a guild.')
    user_id = ctx.author.id
    async with ctx.bot.db_pool.acquire() as conn:
        user_row = await conn.fetchrow(
            'SELECT coordinator_channel_ids,coordinator_room_names FROM users WHERE discord_snowflake=$1',
            user_id
        )
    if not user_row:
        raise NotAValidDatabaseQuery('User not found in database.')
    c_chan = user_row.get('coordinator_channel_ids') or []
    c_room = user_row.get('coordinator_room_names') or []
    room_name = getattr(ctx.channel, 'name', None)
    if target_channel_id not in c_chan and room_name not in c_room:
        raise NotCoordinator('You are not a coordinator in the target channel.')
    return True

async def is_coordinator_via_objects(member, channel):
    if member.guild is None:
        raise NotCoordinator('Command must be used in a guild.')
    async with member._state._get_client().db_pool.acquire() as conn:
        user_row = await conn.fetchrow(
            'SELECT coordinator_channel_ids,coordinator_room_names FROM users WHERE discord_snowflake=$1',
            member.id
        )
    if not user_row:
        return False
    c_chan = user_row.get('coordinator_channel_ids') or []
    c_room = user_row.get('coordinator_room_names') or []
    room_name = getattr(channel, 'name', None)
    if channel.id in c_chan or room_name in c_room:
        return True
    return False
    
async def is_moderator_via_objects(member, channel):
    if member.guild is None:
        raise NotModerator('Command must be used in a guild.')
    async with member._state._get_client().db_pool.acquire() as conn:
        user_row = await conn.fetchrow(
            'SELECT moderator_channel_ids,moderator_room_names FROM users WHERE discord_snowflake=$1',
            member.id
        )
    if not user_row:
        return False
    m_chan = user_row.get('moderator_channel_ids') or []
    m_room = user_row.get('moderator_room_names') or []
    room_name = getattr(channel, 'name', None)
    if channel.id in m_chan or room_name in m_room:
        return True
    return False
    
async def is_owner_developer_via_objects(member: discord.Member, bot: commands.Bot) -> bool:
    for check in (is_system_owner_member, is_guild_owner_member, is_developer_member):
        try:
            if check is is_system_owner_member:
                if await check(member, bot):
                    return True
            else:
                if await check(member):
                    return True
        except Exception:
            continue
    return False

async def is_moderator(ctx, target_channel_id: int):
    if ctx.guild is None:
        raise NotInAGuild('Command must be used in a guild.')
    user_id = ctx.author.id
    room_name = getattr(ctx.channel, 'name', None)
    async with ctx.bot.db_pool.acquire() as conn:
        user_row = await conn.fetchrow(
            'SELECT moderator_channel_ids,moderator_room_names FROM users WHERE discord_snowflake=$1',
            user_id
        )
    if not user_row:
        raise NotAValidDatabaseQuery('User not found in database.')
    m_chan = user_row.get('moderator_channel_ids') or []
    m_room = user_row.get('moderator_room_names') or []
    if target_channel_id not in m_chan and room_name not in m_room:
        raise NotModerator('You are not a moderator in the target channel.')
    return True
    
def is_coordinator_in_channel(channel_id: int):
    async def predicate(ctx: commands.Context):
        room_name = getattr(ctx.channel, 'name', None)
        async with ctx.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT coordinator_channel_ids,coordinator_room_names FROM users WHERE discord_snowflake=$1',
                ctx.author.id
            )
        if row:
            c_chan = row.get('coordinator_channel_ids') or []
            c_room = row.get('coordinator_room_names') or []
            if channel_id in c_chan or room_name in c_room:
                return True
        raise commands.CheckFailure('You are not a coordinator in this channel.')
    return commands.check(predicate)

async def check_block(ctx: commands.Context, member: discord.Member, channel: Optional[discord.abc.GuildChannel] = None):
    bot = ctx.bot
    guild = ctx.guild
    channel_id = getattr(channel, 'id', None) or ctx.channel.id
    room_name = getattr(channel, 'name', None)
    role_hierarchy = ['Everyone', 'Moderator', 'Coordinator', 'Developer', 'Owner']
    author_roles = []
    async with bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT coordinator_channel_ids,moderator_channel_ids,coordinator_room_names,moderator_room_names,developer_guild_ids FROM users WHERE discord_snowflake=$1',
            ctx.author.id
        )
    if row:
        c_chan = row.get('coordinator_channel_ids') or []
        m_chan = row.get('moderator_channel_ids') or []
        c_room = row.get('coordinator_room_names') or []
        m_room = row.get('moderator_room_names') or []
        dev_guilds = row.get('developer_guild_ids') or []
        if channel_id in m_chan or room_name in m_room:
            author_roles.append('Moderator')
        if channel_id in c_chan or room_name in c_room:
            author_roles.append('Coordinator')
        if guild and guild.id in dev_guilds:
            author_roles.append('Developer')
    if guild and ctx.author.id == guild.owner_id:
        author_roles.append('Owner')
    if ctx.author.id == int(bot.config['discord_owner_id']):
        author_roles.append('Owner')
    author_highest = max(author_roles, key=lambda r: role_hierarchy.index(r)) if author_roles else 'Everyone'
    async with bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT coordinator_channel_ids,moderator_channel_ids,coordinator_room_names,moderator_room_names,developer_guild_ids FROM users WHERE discord_snowflake=$1',
            member.id
        )
    target_roles = []
    if row:
        c_chan = row.get('coordinator_channel_ids') or []
        m_chan = row.get('moderator_channel_ids') or []
        c_room = row.get('coordinator_room_names') or []
        m_room = row.get('moderator_room_names') or []
        dev_guilds = row.get('developer_guild_ids') or []
        if channel_id in m_chan or room_name in m_room:
            target_roles.append('Moderator')
        if channel_id in c_chan or room_name in c_room:
            target_roles.append('Coordinator')
        if guild and guild.id in dev_guilds:
            target_roles.append('Developer')
    if guild and member.id == guild.owner_id:
        target_roles.append('Owner')
    if member.id == int(bot.config['discord_owner_id']):
        target_roles.append('Owner')
    target_highest = max(target_roles, key=lambda r: role_hierarchy.index(r)) if target_roles else 'Everyone'
    if member.id == ctx.author.id:
        if 'Owner' in author_roles:
            return 'Owner', True
        return target_highest, False
    success = role_hierarchy.index(target_highest) < role_hierarchy.index(author_highest)
    return target_highest, success
    
async def is_owner_developer_coordinator_via_alias(ctx: commands.Context, alias_type: Optional[str] = None) -> bool:
    guild = ctx.guild
    guild_id = guild.id if guild else None
    current_channel = ctx.channel
    current_channel_id = getattr(current_channel, 'id', None)
    current_room_name = getattr(current_channel, 'name', None)
    guild_aliases = getattr(ctx.bot, 'command_aliases', {}).get(guild_id, {}) if guild_id else {}
    target_channel_id = None
    target_room_name = None
    if alias_type and guild_id:
        alias_data_role = guild_aliases.get('role_aliases', {}).get(alias_type, {})
        alias_data_chan = guild_aliases.get('channel_aliases', {}).get(alias_type, {})
        for alias_data in (alias_data_role, alias_data_chan):
            if alias_data and isinstance(alias_data, dict) and 'channel_id' in alias_data:
                target_channel_id = int(alias_data['channel_id'])
                break
        temp_aliases = guild_aliases.get('temp_room_aliases', {}).get(alias_type, {})
        if temp_aliases:
            alias_name = ctx.command.name.lower()
            if alias_name in temp_aliases:
                target_room_name = temp_aliases[alias_name].get('room_name')

        if not target_channel_id and not target_room_name and ctx.command:
            async with ctx.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(
                    'SELECT channel_id FROM command_aliases WHERE guild_id=$1 AND alias_name=$2',
                    guild_id, ctx.command.name.lower()
                )
            if row:
                target_channel_id = int(row['channel_id'])
    for check in (is_system_owner, is_guild_owner, is_developer):
        try:
            if await check(ctx): return True
        except commands.CheckFailure:
            continue
    async with ctx.bot.db_pool.acquire() as conn:
        user_row = await conn.fetchrow(
            'SELECT coordinator_channel_ids,coordinator_room_names FROM users WHERE discord_snowflake=$1',
            ctx.author.id
        )
    if not user_row:
        return False
    c_chan = user_row.get('coordinator_channel_ids') or []
    c_room = user_row.get('coordinator_room_names') or []
    if target_channel_id in c_chan:
        return True
    if target_room_name in c_room:
        return True
    if current_room_name in c_room:
        return True
    return False

def is_server_muter_predicator():
    async def predicate(ctx: commands.Context):
        guild_id = ctx.guild.id
        user_id = ctx.author.id
        async with ctx.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT server_muter_guild_ids
                FROM users
                WHERE discord_snowflake = $1
            ''', user_id)
        if row and row['server_muter_guild_ids'] and guild_id in row['server_muter_guild_ids']:
            return True
        raise commands.CheckFailure('You do not have server mute permissions in this guild.')
    predicate._permission_level = 'Server Muter'
    return commands.check(predicate)

async def is_owner_app(interaction: discord.Interaction) -> bool:
    system_owner_id = int(interaction.client.config['discord_owner_id'])
    if interaction.user.id == system_owner_id:
        return True
    if interaction.guild and interaction.user.id == interaction.guild.owner_id:
        return True
    return False

async def is_developer_app(interaction: discord.Interaction) -> bool:
    if interaction.guild is None:
        return False
    async with interaction.client.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT developer_guild_ids FROM users WHERE discord_snowflake=$1',
            interaction.user.id
        )
    return row and row.get('developer_guild_ids') and interaction.guild.id in row['developer_guild_ids']

async def is_coordinator_app(interaction: discord.Interaction, alias_type: Optional[str] = None) -> bool:
    if not interaction.guild:
        return False
    target_channel_id = None
    target_room_name = None
    if alias_type:
        guild_aliases = getattr(interaction.client, 'command_aliases', {}).get(interaction.guild.id, {})
        role_alias = guild_aliases.get('role_aliases', {}).get(alias_type, {})
        chan_alias = guild_aliases.get('channel_aliases', {}).get(alias_type, {})
        for a in (role_alias, chan_alias):
            if a and isinstance(a, dict) and 'channel_id' in a:
                target_channel_id = a.get('channel_id')
                break
        temp_aliases = guild_aliases.get('temp_room_aliases', {}).get(alias_type, {})
        if temp_aliases:
            alias_name = interaction.command.name.lower() if interaction.command else None
            if alias_name and alias_name in temp_aliases:
                target_room_name = temp_aliases[alias_name].get('room_name')

    async with interaction.client.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT coordinator_channel_ids, coordinator_room_names FROM users WHERE discord_snowflake=$1',
            interaction.user.id
        )
    if not row:
        return False
    coord_ids = row.get('coordinator_channel_ids') or []
    coord_rooms = row.get('coordinator_room_names') or []
    if alias_type is None and (coord_ids or coord_rooms):
        return True
    if (target_channel_id and target_channel_id in coord_ids) or \
       (target_room_name and target_room_name in coord_rooms):
        return True
    return False

async def is_moderator_app(interaction: discord.Interaction, alias_type: Optional[str] = None) -> bool:
    if not interaction.guild:
        return False
    target_channel_id = None
    target_room_name = None
    if alias_type:
        guild_aliases = getattr(interaction.client, 'command_aliases', {}).get(interaction.guild.id, {})
        role_alias = guild_aliases.get('role_aliases', {}).get(alias_type, {})
        chan_alias = guild_aliases.get('channel_aliases', {}).get(alias_type, {})
        for a in (role_alias, chan_alias):
            if a and isinstance(a, dict) and 'channel_id' in a:
                target_channel_id = a.get('channel_id')
                break
        temp_aliases = guild_aliases.get('temp_room_aliases', {}).get(alias_type, {})
        if temp_aliases:
            alias_name = interaction.command.name.lower() if interaction.command else None
            if alias_name and alias_name in temp_aliases:
                target_room_name = temp_aliases[alias_name].get('room_name')

    async with interaction.client.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT moderator_channel_ids, moderator_room_names FROM users WHERE discord_snowflake=$1',
            interaction.user.id
        )
    if not row:
        return False
    mod_ids = row.get('moderator_channel_ids') or []
    mod_rooms = row.get('moderator_room_names') or []
    if alias_type is None and (mod_ids or mod_rooms):
        return True
    if (target_channel_id and (target_channel_id in mod_ids)) or \
       (target_room_name and (target_room_name in mod_rooms)):
        return True
    return False

def is_owner_app_predicator():
    async def predicate(interaction: discord.Interaction) -> bool:
        if await is_owner_app(interaction):
            return True
        raise app_commands.CheckFailure("You are not system owner or guild owner.")
    return app_commands.check(predicate)

# Owner or Developer
def is_owner_developer_app_predicator():
    async def predicate(interaction: discord.Interaction) -> bool:
        if await is_owner_app(interaction) or await is_developer_app(interaction):
            return True
        raise app_commands.CheckFailure("You are not system owner, guild owner, or developer.")
    return app_commands.check(predicate)

# Owner, Developer, or Coordinator
def is_owner_developer_coordinator_app_predicator(alias_type: Optional[str] = None):
    async def predicate(interaction: discord.Interaction) -> bool:
        if await is_owner_app(interaction) or await is_developer_app(interaction) or await is_coordinator_app(interaction, alias_type):
            return True
        raise app_commands.CheckFailure("You are not system owner, developer, or coordinator.")
    return app_commands.check(predicate)

# Owner, Developer, Coordinator, or Moderator
def is_owner_developer_coordinator_moderator_app_predicator(alias_type: Optional[str] = None):
    async def predicate(interaction: discord.Interaction) -> bool:
        if await is_owner_app(interaction) or await is_developer_app(interaction) \
           or await is_coordinator_app(interaction, alias_type) or await is_moderator_app(interaction, alias_type):
            return True
        raise app_commands.CheckFailure("You are not system owner, developer, coordinator, or moderator.")
    return app_commands.check(predicate)

# Owner, Developer, Coordinator, Moderator, or team role
def is_owner_developer_coordinator_moderator_role_app_predicator(alias_type: Optional[str] = None):
    async def predicate(interaction: discord.Interaction) -> bool:
        # First, check base roles
        if await is_owner_app(interaction) or await is_developer_app(interaction) \
           or await is_coordinator_app(interaction, alias_type) or await is_moderator_app(interaction, alias_type):
            return True

        # Check team roles if applicable
        if interaction.guild:
            db_pool = interaction.client.db_pool
            guild_aliases = getattr(interaction.client, 'command_aliases', {}).get(interaction.guild.id, {})
            target_channel_id = None
            if alias_type:
                role_alias = guild_aliases.get('role_aliases', {}).get(alias_type, {})
                chan_alias = guild_aliases.get('channel_aliases', {}).get(alias_type, {})
                for a in (role_alias, chan_alias):
                    if a and isinstance(a, dict) and 'channel_id' in a:
                        target_channel_id = a.get('channel_id')
                        break
            async with db_pool.acquire() as conn:
                role_rows = await conn.fetch('SELECT role_id FROM role_permissions WHERE is_team_member=TRUE')
            role_ids = [r['role_id'] for r in role_rows] if role_rows else []
            user_role_ids = [role.id for role in getattr(interaction.user, 'roles', [])]
            if any(rid in role_ids for rid in user_role_ids):
                return True

        raise app_commands.CheckFailure("You are not system owner, developer, coordinator, moderator, or in an authorized team role.")
    return app_commands.check(predicate)

# Server muter
def is_server_muter_app_predicator():
    async def predicate(interaction: discord.Interaction) -> bool:
        if not interaction.guild:
            raise app_commands.CheckFailure("This command must be used in a guild.")
        async with interaction.client.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT server_muter_guild_ids FROM users WHERE discord_snowflake=$1',
                interaction.user.id
            )
        if row and row['server_muter_guild_ids'] and interaction.guild.id in row['server_muter_guild_ids']:
            return True
        raise app_commands.CheckFailure("You do not have server mute permissions in this guild.")
    return app_commands.check(predicate)

async def check_owner_dev_coord_mod_app(interaction: discord.Interaction, channel: Optional[discord.abc.GuildChannel] = None) -> Tuple[bool, bool]:
    is_owner_or_dev = False
    is_mod_or_coord = False
    guild = interaction.guild
    guild_id = guild.id if guild else None
    user_id = interaction.user.id
    target_channel_id = getattr(channel, 'id', None) or getattr(interaction, '_target_channel_id', None) or getattr(interaction.channel, 'id', None)
    target_room_name = getattr(channel, 'name', None)
    try:
        if await is_system_owner_app(interaction):
            is_owner_or_dev = True
    except app_commands.CheckFailure:
        try:
            if guild and await is_guild_owner_app(interaction):
                is_owner_or_dev = True
        except app_commands.CheckFailure:
            try:
                if guild and await is_developer_app(interaction):
                    is_owner_or_dev = True
            except app_commands.CheckFailure:
                pass
    if guild:
        async with interaction.client.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids,moderator_channel_ids,coordinator_room_names,moderator_room_names '
                'FROM users WHERE discord_snowflake=$1',
                user_id
            )
        if user_row:
            c_chan = user_row.get('coordinator_channel_ids') or []
            m_chan = user_row.get('moderator_channel_ids') or []
            c_room = user_row.get('coordinator_room_names') or []
            m_room = user_row.get('moderator_room_names') or []
            if (target_channel_id and (target_channel_id in c_chan or target_channel_id in m_chan)) or \
               (target_room_name and (target_room_name in c_room or target_room_name in m_room)):
                is_mod_or_coord = True
    return is_owner_or_dev, is_mod_or_coord

async def is_guild_owner_app(interaction: discord.Interaction) -> bool:
    if interaction.guild is None:
        raise NotInAGuild('Command must be used in a guild.')
    if interaction.guild.owner_id != interaction.user.id:
        raise NotGuildOwner()
    return True

async def is_system_owner_app(interaction: discord.Interaction) -> bool:
    system_owner_id = int(interaction.client.config['discord_owner_id'])
    if interaction.user.id != system_owner_id:
        raise NotSystemOwner()
    return True

async def check_block_app(ctx: discord.Interaction, member: discord.Member, channel: Optional[discord.abc.GuildChannel]=None):
    bot=ctx.client
    channel_id=getattr(channel,'id',None) or ctx.channel.id
    guild=ctx.guild
    role_hierarchy=['Everyone','Moderator','Coordinator','Developer','Owner']
    room_name=getattr(channel,'name',None)
    async with bot.db_pool.acquire() as conn:
        row=await conn.fetchrow('SELECT coordinator_channel_ids,moderator_channel_ids,coordinator_room_names,moderator_room_names,developer_guild_ids FROM users WHERE discord_snowflake=$1',ctx.user.id)
    author_roles=[]
    if row:
        c_chan=row.get('coordinator_channel_ids') or []
        c_room=row.get('coordinator_room_names') or []
        dev_guilds=row.get('developer_guild_ids') or []
        m_chan=row.get('moderator_channel_ids') or []
        m_room=row.get('moderator_room_names') or []
        if channel_id in m_chan or room_name in m_room: author_roles.append('Moderator')
        if channel_id in c_chan or room_name in c_room: author_roles.append('Coordinator')
        if guild and guild.id in dev_guilds: author_roles.append('Developer')
    if guild and ctx.user.id==guild.owner_id: author_roles.append('Owner')
    if ctx.user.id==int(bot.config['discord_owner_id']): author_roles.append('Owner')
    author_highest=max(author_roles,key=lambda r: role_hierarchy.index(r)) if author_roles else 'Everyone'
    async with bot.db_pool.acquire() as conn:
        row=await conn.fetchrow('SELECT coordinator_channel_ids,moderator_channel_ids,coordinator_room_names,moderator_room_names,developer_guild_ids FROM users WHERE discord_snowflake=$1',member.id)
    target_roles=[]
    if row:
        c_chan=row.get('coordinator_channel_ids') or []
        c_room=row.get('coordinator_room_names') or []
        dev_guilds=row.get('developer_guild_ids') or []
        m_chan=row.get('moderator_channel_ids') or []
        m_room=row.get('moderator_room_names') or []
        if channel_id in m_chan or room_name in m_room: target_roles.append('Moderator')
        if channel_id in c_chan or room_name in c_room: target_roles.append('Coordinator')
        if guild and guild.id in dev_guilds: target_roles.append('Developer')
    if guild and member.id==guild.owner_id: target_roles.append('Owner')
    if member.id==int(bot.config['discord_owner_id']): target_roles.append('Owner')
    target_highest=max(target_roles,key=lambda r: role_hierarchy.index(r)) if target_roles else 'Everyone'
    if member.id==ctx.user.id:
        if 'Owner' in author_roles: return 'Owner',True
        return target_highest,False
    success=role_hierarchy.index(target_highest)<role_hierarchy.index(author_highest)
    return target_highest,success

async def check_owner_dev_coord_mod_overall(ctx: commands.Context) -> Tuple[bool,bool]:
    guild = ctx.guild
    is_mod_or_coord = False
    is_owner_or_dev = False
    user_id = ctx.author.id
    try:
        if await is_system_owner(ctx): is_owner_or_dev = True
    except commands.CheckFailure:
        try:
            if guild and await is_guild_owner(ctx): is_owner_or_dev = True
        except commands.CheckFailure:
            try:
                if guild and await is_developer(ctx): is_owner_or_dev = True
            except commands.CheckFailure: pass
    if guild:
        async with ctx.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT coordinator_channel_ids,coordinator_room_names,moderator_channel_ids,moderator_room_names FROM users WHERE discord_snowflake=$1',user_id)
        if row:
            if (row.get('coordinator_channel_ids') or []) or (row.get('coordinator_room_names') or []) or (row.get('moderator_channel_ids') or []) or (row.get('moderator_room_names') or []):
                is_mod_or_coord = True
    return is_owner_or_dev,is_mod_or_coord

async def check_owner_dev_coord_overall(ctx: commands.Context) -> Tuple[bool,bool]:
    guild = ctx.guild
    is_coord = False
    is_owner_or_dev = False
    user_id = ctx.author.id
    try:
        if await is_system_owner(ctx): is_owner_or_dev = True
    except commands.CheckFailure:
        try:
            if guild and await is_guild_owner(ctx): is_owner_or_dev = True
        except commands.CheckFailure:
            try:
                if guild and await is_developer(ctx): is_owner_or_dev = True
            except commands.CheckFailure: pass
    if guild:
        async with ctx.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT coordinator_channel_ids,coordinator_room_names FROM users WHERE discord_snowflake=$1',user_id)
        if row:
            if (row.get('coordinator_channel_ids') or []) or (row.get('coordinator_room_names') or []):
                is_coord = True
    return is_owner_or_dev,is_coord

async def check_owner_dev_coord_mod_overall_app(interaction: discord.Interaction) -> Tuple[bool,bool]:
    guild = interaction.guild
    is_mod_or_coord = False
    is_owner_or_dev = False
    user_id = interaction.user.id
    try:
        if await is_system_owner_inter(interaction): is_owner_or_dev = True
    except app_commands.CheckFailure:
        try:
            if guild and await is_guild_owner_inter(interaction): is_owner_or_dev = True
        except app_commands.CheckFailure:
            try:
                if guild and await is_developer_inter(interaction): is_owner_or_dev = True
            except app_commands.CheckFailure: pass
    if guild:
        async with interaction.client.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT coordinator_channel_ids,coordinator_room_names,moderator_channel_ids,moderator_room_names FROM users WHERE discord_snowflake=$1',user_id)
        if row:
            if (row.get('coordinator_channel_ids') or []) or (row.get('coordinator_room_names') or []) or (row.get('moderator_channel_ids') or []) or (row.get('moderator_room_names') or []):
                is_mod_or_coord = True
    return is_owner_or_dev,is_mod_or_coord

async def check_owner_dev_coord_overall_app(interaction: discord.Interaction) -> Tuple[bool,bool]:
    guild = interaction.guild
    is_coord = False
    is_owner_or_dev = False
    user_id = interaction.user.id
    try:
        if await is_system_owner_inter(interaction): is_owner_or_dev = True
    except app_commands.CheckFailure:
        try:
            if guild and await is_guild_owner_inter(interaction): is_owner_or_dev = True
        except app_commands.CheckFailure:
            try:
                if guild and await is_developer_inter(interaction): is_owner_or_dev = True
            except app_commands.CheckFailure: pass
    if guild:
        async with interaction.client.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT coordinator_channel_ids,coordinator_room_names FROM users WHERE discord_snowflake=$1',user_id)
        if row:
            if (row.get('coordinator_channel_ids') or []) or (row.get('coordinator_room_names') or []):
                is_coord = True
    return is_owner_or_dev,is_coord
