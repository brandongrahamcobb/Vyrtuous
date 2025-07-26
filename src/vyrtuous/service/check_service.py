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

import discord
from discord.ext import commands
from vyrtuous.utils.setup_logging import logger as permission_logger
from vyrtuous.config import Config  # adjust import path accordingly

config = Config.get_config()
      
class NotCoordinator(commands.CheckFailure):
    def __init__(self, message="You are not a coordinator in the requested channel."):
        super().__init__(message)

class NotModerator(commands.CheckFailure):
    def __init__(self, message="You are not a moderator in the requested channel."):
        super().__init__(message)

class NotGuildOwner(commands.CheckFailure):
    def __init__(self, message="You are not the guild owner."):
        super().__init__(message)

class NotSystemOwner(commands.CheckFailure):
    def __init__(self, message="You are not the bot owner."):
        super().__init__(message)

class NotDeveloper(commands.CheckFailure):
    def __init__(self, message="You are not a developer in this guild."):
        super().__init__(message)

class NoCommandAlias(commands.CheckFailure):
    def __init__(self, message="This command alias is not mapped to a target channel."):
        super().__init__(message)
                                    
class NotAtHome(commands.CheckFailure):
    def __init__(self, message="You are not in the home Discord!."):
        super().__init__(message)

async def at_home(ctx) -> bool:
    if ctx.guild is not None and ctx.guild.id == config['discord_testing_guild_id']:
        return True
    raise NotAtHome()
                                    
async def has_command_alias(ctx) -> bool:
    if not ctx.guild or not ctx.command:
        raise NoCommandAlias("Command must be used in a guild and must be a valid command.")
    bot = ctx.bot
    guild_id = ctx.guild.id
    command_name = ctx.command.name.lower()       # e.g. 'vmute'
    invoked_name = ctx.invoked_with.lower()       # what the user actually typed
    for alias_type in ("mute", "unmute", "ban", "unban"):
        alias_map = bot.command_aliases.get(guild_id, {}).get(alias_type, {})
        if command_name in alias_map:
            target_channel_id = alias_map[command_name]
            ctx._target_channel_id = target_channel_id
            return True
    raise NoCommandAlias(f"Command alias `{command_name}` is not mapped to any target channel.")

async def is_moderator(ctx):
    user_id = ctx.author.id
    guild_id = ctx.guild.id if ctx.guild else None
    channel_id = ctx.channel.id
    
    permission_logger.debug(f"=== MODERATOR CHECK START ===")
    permission_logger.debug(f"User ID: {user_id}")
    permission_logger.debug(f"Guild ID: {guild_id}")
    permission_logger.debug(f"Channel ID: {channel_id}")
    
    if ctx.guild is None:
        permission_logger.debug("❌ No guild context")
        raise commands.CheckFailure("Command must be used in a guild.")
    
    bot = ctx.bot
    
    try:
        async with bot.db_pool.acquire() as conn:
            permission_logger.debug("📊 Executing database query...")
            row = await conn.fetchrow(
                "SELECT moderator_ids, moderator_channel_ids FROM users WHERE user_id = $1",
                user_id
            )
            permission_logger.debug(f"📊 Database row result: {row}")
            
    except Exception as e:
        permission_logger.error(f"❌ Database error: {e}")
        raise commands.CheckFailure(f"Database error: {str(e)}")
    
    if not row:
        permission_logger.debug("❌ User not found in database")
        raise commands.CheckFailure("User not found in database.")
    
    moderator_ids = row.get("moderator_ids") if row else None
    moderator_channels = row.get("moderator_channel_ids") if row else None
    
    permission_logger.debug(f"🔑 Moderator IDs from DB: {moderator_ids}")
    permission_logger.debug(f"📺 Moderator channels from DB: {moderator_channels}")
    
async def is_coordinator(ctx):
    user_id = ctx.author.id
    guild_id = ctx.guild.id if ctx.guild else None
    channel_id = ctx.channel.id
    
    permission_logger.debug(f"=== COORDINATOR CHECK START ===")
    permission_logger.debug(f"User ID: {user_id}")
    permission_logger.debug(f"Guild ID: {guild_id}")
    permission_logger.debug(f"Channel ID: {channel_id}")
    
    if ctx.guild is None:
        permission_logger.debug("❌ No guild context")
        raise commands.CheckFailure("Command must be used in a guild.")
    
    bot = ctx.bot
    
    try:
        async with bot.db_pool.acquire() as conn:
            permission_logger.debug("📊 Executing database query...")
            row = await conn.fetchrow(
                "SELECT coordinator_ids, coordinator_channel_ids FROM users WHERE user_id = $1",
                user_id
            )
            permission_logger.debug(f"📊 Database row result: {row}")
            
    except Exception as e:
        permission_logger.error(f"❌ Database error: {e}")
        raise commands.CheckFailure(f"Database error: {str(e)}")
    
    if not row:
        permission_logger.debug("❌ User not found in database")
        raise commands.CheckFailure("User not found in database.")
    
    coordinator_ids = row.get("coordinator_ids") if row else None
    coordinator_channels = row.get("coordinator_channel_ids") if row else None
    
    permission_logger.debug(f"🔑 Coordinator IDs from DB: {coordinator_ids}")
    permission_logger.debug(f"📺 Coordinator channels from DB: {coordinator_channels}")
    
    # Check guild permissions
    if not coordinator_ids:
        permission_logger.debug("❌ No coordinator_ids found for user")
        raise commands.CheckFailure("You have no coordinator permissions configured.")
    
    if guild_id not in coordinator_ids:
        permission_logger.debug(f"❌ Guild ID {guild_id} not in coordinator_ids {coordinator_ids}")
        raise commands.CheckFailure("You are not a coordinator in this guild.")
    
    permission_logger.debug(f"✅ Guild permission check passed")
    
    # Check channel permissions
    if not coordinator_channels:
        permission_logger.debug("❌ No coordinator_channel_ids found for user")
        raise commands.CheckFailure("You have no channel permissions configured.")
    
    if channel_id not in coordinator_channels:
        permission_logger.debug(f"❌ Channel ID {channel_id} not in coordinator_channel_ids {coordinator_channels}")
        raise commands.CheckFailure("You are not authorized to use commands in this channel.")
    
    permission_logger.debug(f"✅ Channel permission check passed")
    permission_logger.debug("=== COORDINATOR CHECK SUCCESS ===")
    return True
async def is_developer(ctx):
    if ctx.guild is None:
        raise NotDeveloper("Command must be used in a guild.")
    async with ctx.bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT developer_guild_ids FROM users WHERE user_id = $1", ctx.author.id
        )
    if not row or not row.get("developer_guild_ids") or ctx.guild.id not in row.get("developer_guild_ids", []):
        raise NotDeveloper("You are not a developer in this guild.")
    return True
                                    
async def is_guild_owner(ctx):
    if ctx.guild is None:
        raise NotGuildOwner("Command must be used in a guild.")
    if ctx.guild.owner_id != ctx.author.id:
        raise NotGuildOwner()
    return True
                    
async def is_owner(ctx):
    errors = []
    for check in (is_guild_owner, is_system_owner):
        try:
            if await check(ctx):
                return True
        except commands.CheckFailure as e:
            errors.append(str(e))
    raise commands.CheckFailure("\n".join(f"❌ {msg}" for msg in errors))

async def is_owner_developer(ctx):
    errors = []
    for check in (is_guild_owner, is_system_owner, is_developer):
        try:
            if await check(ctx):
                return True
        except commands.CheckFailure as e:
            errors.append(str(e))
    raise commands.CheckFailure("\n".join(f"❌ {msg}" for msg in errors))
    
async def is_owner_developer_coordinator(ctx):
    errors = []
    for check in (is_guild_owner, is_system_owner, is_developer, is_coordinator):
        try:
            if await check(ctx):
                return True
        except commands.CheckFailure as e:
            errors.append(str(e))
    raise commands.CheckFailure("\n".join(f"❌ {msg}" for msg in errors))
    
async def is_owner_developer_coordinator_moderator(ctx):
    errors = []
    for check in (is_guild_owner, is_system_owner, is_developer, is_coordinator, is_moderator):
        try:
            if await check(ctx):
                return True
        except commands.CheckFailure as e:
            errors.append(str(e))
    raise commands.CheckFailure("\n".join(f"❌ {msg}" for msg in errors))
    
async def is_system_owner(ctx):
    system_owner_id = ctx.bot.config["discord_owner_id"]
    if ctx.author.id != system_owner_id:
        raise NotSystemOwner()
    return True
                                    

