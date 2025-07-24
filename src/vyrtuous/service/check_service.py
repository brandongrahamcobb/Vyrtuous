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
from vyrtuous.utils.setup_logging import logger
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
    for alias_type in ("mute", "unmute"):
        alias_map = bot.command_aliases.get(guild_id, {}).get(alias_type, {})
        if command_name in alias_map:
            target_channel_id = alias_map[command_name]
            ctx._target_channel_id = target_channel_id
            return True
    raise NoCommandAlias(f"Command alias `{command_name}` is not mapped to any target channel.")

async def is_channel_moderator(ctx):
    bot = ctx.bot
    user_id = ctx.author.id
    guild_id = ctx.guild.id
    command_name = ctx.invoked_with.lower()
    target_channel_id = None
    for alias_type in ("mute", "unmute"):
        alias_map = bot.command_aliases.get(guild_id, {}).get(alias_type, {})
        if command_name in alias_map:
            target_channel_id = alias_map.get(command_name)
            break
    if not target_channel_id:
        raise NotModerator("This command alias is not mapped to a target channel.")
    async with bot.db_pool.acquire() as conn:
        row = await conn.fetchrow("""
            SELECT moderator_ids FROM users WHERE user_id = $1
        """, user_id)
    if not row:
        raise NotModerator(f"You are not a VC moderator in <#{target_channel_id}>.")
    moderator_ids = row.get("moderator_ids") or []
    if target_channel_id not in moderator_ids:
        raise NotModerator()
    return True

async def is_coordinator(ctx):
    if ctx.guild is None:
        raise NotCoordinator("Command must be used in a guild.")
    async with ctx.bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT coordinator_ids FROM users WHERE user_id = $1", ctx.author.id
        )
    if not row or ctx.guild.id not in row.get("coordinator_ids", []):
        raise NotCoordinator()
    return True
    
async def is_developer(ctx):
    if ctx.guild is None:
        raise NotDeveloper("Command must be used in a guild.")
    async with ctx.bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT developer_guild_ids FROM users WHERE user_id = $1", ctx.author.id
        )
    if not row or ctx.guild.id not in row.get("developer_guild_ids", []):
        raise NotDeveloper()
    return True
                                    
async def is_guild_owner(ctx):
    if ctx.guild is None:
        raise NotGuildOwner("Command must be used in a guild.")
    if ctx.guild.owner_id != ctx.author.id:
        raise NotGuildOwner()
    return True
                                    
async def is_moderator(ctx):
    errors = []
    for check in (has_command_alias, is_channel_moderator):
        try:
            if await check(ctx):
                continue
        except commands.CheckFailure as e:
            errors.append(str(e))
    if errors:
        raise commands.CheckFailure("\n".join(f"❌ {e}" for e in errors))
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
                                    

