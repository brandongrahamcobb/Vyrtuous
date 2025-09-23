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
    for alias_type in ('role', 'unrole', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute', 'cow', 'uncow'):
        alias_map = bot.command_aliases.get(guild_id, {}).get(alias_type, {})
        if command_name in alias_map:
            target_channel_id = alias_map[command_name]
            ctx._target_channel_id = target_channel_id
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
        target_channel_id = None
        current_channel_id = getattr(ctx.channel, 'id', None)
        if alias_type and ctx.guild:
            guild_aliases = getattr(ctx.bot, 'command_aliases', {}).get(ctx.guild.id, {})
            alias_data_role = guild_aliases.get('role_aliases', {}).get(alias_type, {})
            alias_data_chan = guild_aliases.get('channel_aliases', {}).get(alias_type, {})
            for alias_data in (alias_data_role, alias_data_chan):
                if alias_data and alias_data.get('channel_id') == current_channel_id:
                    target_channel_id = alias_data['channel_id']
                    break
        if not target_channel_id and ctx.guild and ctx.command:
            async with ctx.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(
                    'SELECT channel_id FROM command_aliases WHERE guild_id = $1 AND alias_name = $2',
                    ctx.guild.id,
                    ctx.command.name.lower()
                )
            if row and row['channel_id'] == current_channel_id:
                target_channel_id = row['channel_id']
        if not target_channel_id:
            return True
        for check in (is_system_owner, is_guild_owner, is_developer):
            try:
                if await check(ctx): return True
            except commands.CheckFailure: continue
        async with ctx.bot.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids, moderator_channel_ids FROM users WHERE discord_snowflake = $1',
                ctx.author.id
            )
        if user_row:
            if target_channel_id in (user_row.get('coordinator_channel_ids') or []): return True
            if target_channel_id in (user_row.get('moderator_channel_ids') or []): return True
        raise commands.CheckFailure('You are not system owner, guild owner, developer' + (', coordinator' if alias_type else '') + (', or moderator' if alias_type else ''))
    predicate._permission_level = 'Moderator'
    return commands.check(predicate)


def is_owner_developer_coordinator_predicator(alias_type: Optional[str] = None):
    async def predicate(ctx: commands.Context):
        target_channel_id = None
        current_channel_id = getattr(ctx.channel, 'id', None)
        if alias_type and ctx.guild:
            guild_aliases = getattr(ctx.bot, 'command_aliases', {}).get(ctx.guild.id, {})
            alias_data_role = guild_aliases.get('role_aliases', {}).get(alias_type, {})
            alias_data_chan = guild_aliases.get('channel_aliases', {}).get(alias_type, {})
            for alias_data in (alias_data_role, alias_data_chan):
                if alias_data and alias_data.get('channel_id') == current_channel_id:
                    target_channel_id = alias_data['channel_id']
                    break
        if not target_channel_id and ctx.guild and ctx.command:
            async with ctx.bot.db_pool.acquire() as conn:
                row = await conn.fetchrow(
                    'SELECT channel_id FROM command_aliases WHERE guild_id = $1 AND alias_name = $2',
                    ctx.guild.id,
                    ctx.command.name.lower()
                )
            if row and row['channel_id'] == current_channel_id:
                target_channel_id = row['channel_id']
        if not target_channel_id:
            return True
        for check in (is_system_owner, is_guild_owner, is_developer):
            try:
                if await check(ctx): return True
            except commands.CheckFailure: continue
        async with ctx.bot.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1',
                ctx.author.id
            )
        if user_row and target_channel_id in (user_row.get('coordinator_channel_ids') or []):
            return True
        raise commands.CheckFailure('You are not system owner, guild owner, developer, or coordinator in this channel')
    predicate._permission_level = 'Coordinator'
    return commands.check(predicate)

async def check_owner_dev_coord_mod(ctx: commands.Context, channel: Optional[discord.abc.GuildChannel] = None) -> Tuple[bool, bool]:
    is_owner_or_dev = False
    is_mod_or_coord = False
    guild_id = ctx.guild.id if ctx.guild else None
    user_id = ctx.author.id
    target_channel_id = getattr(channel, 'id', None) or getattr(ctx, '_target_channel_id', None) or getattr(ctx.channel, 'id', None)
    try:
        if await is_system_owner(ctx):
            is_owner_or_dev = True
    except commands.CheckFailure:
        if ctx.guild:
            try:
                if await is_guild_owner(ctx):
                    is_owner_or_dev = True
            except commands.CheckFailure:
                if ctx.guild:
                    try:
                        if await is_developer(ctx):
                            is_owner_or_dev = True
                    except commands.CheckFailure:
                        pass
    if ctx.guild and target_channel_id:
        async with ctx.bot.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids, moderator_channel_ids FROM users WHERE discord_snowflake = $1',
                user_id
            )
        if user_row:
            coordinator_ids = user_row.get('coordinator_channel_ids') or []
            moderator_ids = user_row.get('moderator_channel_ids') or []
            if target_channel_id in coordinator_ids:
                is_mod_or_coord = True
            elif target_channel_id in moderator_ids:
                is_mod_or_coord = True

    return is_owner_or_dev, is_mod_or_coord

async def check_owner_dev_coord(ctx: commands.Context, channel: Optional[discord.abc.GuildChannel] = None) -> Tuple[bool, bool]:
    is_owner_or_dev = False
    is_coord = False
    user_id = ctx.author.id
    target_channel_id = getattr(channel, 'id', None) or getattr(ctx, '_target_channel_id', None) or getattr(ctx.channel, 'id', None)
    try:
        if await is_system_owner(ctx):
            is_owner_or_dev = True
    except commands.CheckFailure:
        try:
            if ctx.guild and await is_guild_owner(ctx):
                is_owner_or_dev = True
        except commands.CheckFailure:
            try:
                if ctx.guild and await is_developer(ctx):
                    is_owner_or_dev = True
            except commands.CheckFailure:
                pass
    if ctx.guild and target_channel_id:
        async with ctx.bot.db_pool.acquire() as conn:
            user_row = await conn.fetchrow(
                'SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1',
                user_id
            )
        if user_row:
            coordinator_ids = user_row.get('coordinator_channel_ids') or []
            if target_channel_id in coordinator_ids:
                is_coord = True

    return is_owner_or_dev, is_coord

async def is_coordinator(ctx, target_channel_id: int):
    if ctx.guild is None:
        raise NotCoordinator('Command must be used in a guild.')
    user_id = ctx.author.id
    async with ctx.bot.db_pool.acquire() as conn:
        user_row = await conn.fetchrow(
            'SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1',
            user_id
        )
    if not user_row:
        raise NotAValidDatabaseQuery('User not found in database.')
    coordinator_channel_ids = user_row.get('coordinator_channel_ids') or []
    if target_channel_id not in coordinator_channel_ids:
        raise NotCoordinator('You are not a coordinator in the target channel.')
    return True

async def is_moderator(ctx, target_channel_id: int):
    if ctx.guild is None:
        raise NotInAGuild('Command must be used in a guild.')
    user_id = ctx.author.id
    async with ctx.bot.db_pool.acquire() as conn:
        user_row = await conn.fetchrow(
            'SELECT moderator_channel_ids FROM users WHERE discord_snowflake = $1',
            user_id
        )
    if not user_row:
        raise NotAValidDatabaseQuery('User not found in database.')
    moderator_channel_ids = user_row.get('moderator_channel_ids') or []
    if target_channel_id not in moderator_channel_ids:
        raise NotModerator('You are not a moderator in the target channel.')
    return True
    
def is_coordinator_in_channel(channel_id: int):
    async def predicate(ctx: commands.Context):
        async with ctx.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1',
                ctx.author.id
            )
        if row and channel_id in (row.get('coordinator_channel_ids') or []):
            return True
        raise commands.CheckFailure('You are not a coordinator in this channel.')
    return commands.check(predicate)

async def check_block(ctx: commands.Context, member: discord.Member, channel_id: int):
    bot = ctx.bot
    role_hierarchy = ['Everyone', 'Moderator', 'Coordinator', 'Developer', 'Owner']
    author_roles = []
    async with bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT coordinator_channel_ids, moderator_channel_ids, developer_guild_ids FROM users WHERE discord_snowflake = $1',
            ctx.author.id
        )
    if row:
        if row.get('moderator_channel_ids') and channel_id in row['moderator_channel_ids']:
            author_roles.append('Moderator')
        if row.get('coordinator_channel_ids') and channel_id in row['coordinator_channel_ids']:
            author_roles.append('Coordinator')
        if row.get('developer_guild_ids') and ctx.guild and ctx.guild.id in row['developer_guild_ids']:
            author_roles.append('Developer')
    if ctx.guild:
        if ctx.author.id == ctx.guild.owner_id:
            author_roles.append('Owner')
    if ctx.author.id == int(bot.config['discord_owner_id']):
        author_roles.append('Owner')
    author_highest = max(author_roles, key=lambda r: role_hierarchy.index(r)) if author_roles else 'Everyone'
    target_roles = []
    async with bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            'SELECT coordinator_channel_ids, moderator_channel_ids, developer_guild_ids FROM users WHERE discord_snowflake = $1',
            member.id
        )
    if row:
        if row.get('moderator_channel_ids') and channel_id in row['moderator_channel_ids']:
            target_roles.append('Moderator')
        if row.get('coordinator_channel_ids') and channel_id in row['coordinator_channel_ids']:
            target_roles.append('Coordinator')
        if row.get('developer_guild_ids') and ctx.guild and ctx.guild.id in row['developer_guild_ids']:
            target_roles.append('Developer')
    if ctx.guild:
        if member.id == ctx.author.id:
            if 'Owner' in author_roles:
                return 'Owner', True
            return target_highest, False
    if member.id == int(bot.config['discord_owner_id']):
        target_roles.append('Owner')
    target_highest = max(target_roles, key=lambda r: role_hierarchy.index(r)) if target_roles else 'Everyone'
    success = role_hierarchy.index(target_highest) < role_hierarchy.index(author_highest)
    return target_highest, success
    
async def is_owner_developer_coordinator_via_alias(ctx: commands.Context, alias_type: Optional[str] = None) -> bool:
    print("test")
    current_channel_id = getattr(ctx.channel, 'id', None)
    guild_id = getattr(getattr(ctx, 'guild', None), 'id', None)
    guild_aliases = getattr(ctx.bot, 'command_aliases', {}).get(guild_id, {}) if guild_id else {}
    alias_data_role = guild_aliases.get('role_aliases', {}).get(alias_type, {}) or {}
    alias_data_chan = guild_aliases.get('channel_aliases', {}).get(alias_type, {}) or {}
    target_channel_id = None
    for source, is_role in ((alias_data_role, True), (alias_data_chan, False)):
        for alias_name, alias_data in source.items():
            channel_id = alias_data.get('channel_id') if is_role and isinstance(alias_data, dict) else alias_data
            if channel_id:
                target_channel_id = int(channel_id)
                break
        if target_channel_id: break
    if not target_channel_id and guild_id and ctx.command:
        async with ctx.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT channel_id FROM command_aliases WHERE guild_id = $1 AND alias_name = $2',
                guild_id, ctx.command.name.lower()
            )
        if row:
            target_channel_id = int(row['channel_id'])
    for check in (is_system_owner, is_guild_owner, is_developer):
        try:
            if await check(ctx):
                return True
            else:
                logger.info(f"Permission check {check.__name__} returned False")
        except commands.CheckFailure:
            logger.info(f"Permission check {check.__name__} raised CheckFailure")
            continue
    async with ctx.bot.db_pool.acquire() as conn:
        user_row = await conn.fetchrow(
            'SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1',
            ctx.author.id
        )
    user_channels = user_row['coordinator_channel_ids'] if user_row else []
    if target_channel_id and target_channel_id in (user_channels or []):
        return True
    return False
