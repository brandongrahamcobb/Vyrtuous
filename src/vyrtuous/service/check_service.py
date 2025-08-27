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
from vyrtuous.config import Config
from vyrtuous.bot.discord_bot import DiscordBot# adjust import path accordingly
import discord
config = Config.get_config()
      
class NotCoordinator(commands.CheckFailure):
    def __init__(self, message="You are not a coordinator in the requested channel."):
        super().__init__(message)

class NoInheritedPermissionsForRole(commands.CheckFailure):
    def __init__(self, message="You are not a coordinator in the requested channel."):
        super().__init__(message)

class NotModerator(commands.CheckFailure):
    def __init__(self, message="You are not a moderator in the requested channel."):
        super().__init__(message)

class NotGuildOwner(commands.CheckFailure):
    def __init__(self, message="You are not the guild owner and cannot do this."):
        super().__init__(message)

class CantChoseGuildOwner(commands.CheckFailure):
    def __init__(self, message="You selected the guild owner and cannot do this."):
        super().__init__(message)

class NotInAGuild(commands.CheckFailure):
    def __init__(self, message="You did not send your message in a guild and cannot do this."):
        super().__init__(message)

class NotAValidDatabaseQuery(commands.CheckFailure):
    def __init__(self, message="Database query is invalid."):
        super().__init__(message)

class NoPermissionsForRoleSetup(commands.CheckFailure):
    def __init__(self, message="No permissions for this role are setup."):
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
    command_name = ctx.command.name.lower()
    for alias_type in ("mute", "unmute", "ban", "unban", "flag", "unflag", "tmute", "untmute"):
        alias_map = bot.command_aliases.get(guild_id, {}).get(alias_type, {})
        if command_name in alias_map:
            target_channel_id = alias_map[command_name]
            ctx._target_channel_id = target_channel_id
            return True
    raise NoCommandAlias(f"Command alias `{command_name}` is not mapped to any target channel.")

async def is_moderator(ctx):
    if ctx.guild is None:
        raise NotInAGuild("Command must be used in a guild.")
    user_id = ctx.author.id
    guild_id = ctx.guild.id
    async with ctx.bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT moderator_ids FROM users WHERE user_id = $1", user_id
        )
    if not row or guild_id not in (row.get("moderator_ids") or []):
        raise NotModerator("You are not a moderator in this guild.")
    return True

    # if not row:
    #     raise NotAValidDatabaseQuery("User not found in database.")
    # moderator_ids = row.get("moderator_ids")
    # moderator_channels = row.get("moderator_channel_ids")
    # if not moderator_ids:
    #     raise NoPermissionsForRoleSetup("There are no moderator permissions configured.")
    # if guild_id not in moderator_ids:
    #     raise NotModerator("You are not a moderator in this guild.")
    # if not moderator_channels:
    #     raise NoInheritedPermissionsForRole("You have no channel permissions configured.")
    # return True

async def is_coordinator(ctx):
    if ctx.guild is None:
        raise NotCoordinator("Command must be used in a guild.")
    user_id = ctx.author.id
    guild_id = ctx.guild.id
    bot = ctx.bot
    try:
        async with bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                "SELECT coordinator_ids FROM users WHERE user_id = $1",
                user_id
            )
    except Exception:
        raise NotAValidDatabaseQuery("Failed to query coordinator permissions.")
    if not row:
        raise NotAValidDatabaseQuery("User not found in database.")
    coordinator_ids = row.get("coordinator_ids")
    if not coordinator_ids or guild_id not in coordinator_ids:
        raise NotCoordinator("You are not a coordinator in this guild.")
    return True
# async def is_coordinator(ctx):
#     if ctx.guild is None:
#         raise NotCoordinator("Command must be used in a guild.")
#     user_id = ctx.author.id
#     guild_id = ctx.guild.id
#     channel_id = ctx.channel.id
#     bot = ctx.bot
#     try:
#         async with bot.db_pool.acquire() as conn:
#             row = await conn.fetchrow(
#                 "SELECT coordinator_ids, coordinator_channel_ids FROM users WHERE user_id = $1",
#                 user_id
#             )
#     except Exception as e:
#         raise NotAValidDatabaseQuery()
#     if not row:
#         raise NotAValidDatabaseQuery()
#     coordinator_ids = row.get("coordinator_ids")
#     coordinator_channels = row.get("coordinator_channel_ids")
#     if not coordinator_ids:
#         raise NoPermissionsForRoleSetup("You have no coordinator permissions configured.")
#     if guild_id not in coordinator_ids:
#         raise NotCoordinator("You are not a coordinator in this guild.")
#     if not coordinator_channels:
#         raise NoInheritedPermissionsForRole("You have no channel permissions configured.")
#     if channel_id not in coordinator_channels:
#         raise NotCoordinator("You are not authorized to use coordinator commands in this channel or the specified channel.")
#     return True

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
        raise NotInAGuild("Command must be used in a guild.")
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
    raise commands.CheckFailure("\n".join(f"🔥 {msg}" for msg in errors))

async def is_coordinator_or_moderator_for_channel(ctx: commands.Context, channel: discord.abc.GuildChannel) -> bool:
    user_id = ctx.author.id
    guild_id = ctx.guild.id
    channel_id = channel.id
    async with ctx.bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT coordinator_channel_ids, moderator_ids, moderator_channel_ids FROM users WHERE user_id = $1",
            user_id
        )
    if not row:
        raise NotAValidDatabaseQuery("User not found in database.")
    coordinator_channel_ids = set(row.get("coordinator_channel_ids") or [])
    moderator_channel_ids = set(row.get("moderator_channel_ids") or [])
    moderator_ids = set(row.get("moderator_ids") or [])
    if channel_id in coordinator_channel_ids:
        return True
    if guild_id in moderator_ids and channel_id in moderator_channel_ids:
        return True
    raise NoInheritedPermissionsForRole("You are not a coordinator or moderator for this channel.")

async def is_coordinator_for_channel(ctx: commands.Context, channel: discord.VoiceChannel) -> bool:
    async with ctx.bot.db_pool.acquire() as conn:
        row = await conn.fetchrow('''
            SELECT 1
            FROM users
            WHERE user_id = $1
              AND $2 = ANY(coordinator_channel_ids)
        ''', ctx.author.id, channel.id)
    return row is not None

async def is_moderator_for_channel(ctx, channel: discord.abc.GuildChannel):
    user_id = ctx.author.id
    guild_id = ctx.guild.id
    channel_id = channel.id
    async with ctx.bot.db_pool.acquire() as conn:
        row = await conn.fetchrow(
            "SELECT moderator_ids, moderator_channel_ids FROM users WHERE user_id = $1",
            user_id
        )
    if not row:
        raise NotAValidDatabaseQuery("User not found in database.")
    moderator_ids = set(row.get("moderator_ids") or [])
    channel_ids = set(row.get("moderator_channel_ids") or [])
    if guild_id not in moderator_ids:
        raise NotModerator("You are not a moderator in this guild.")
    if channel_id not in channel_ids:
        raise NotModerator("You are not a moderator in this channel.")
    return True

async def is_guild_owner_block(ctx, user_id: int):
    if user_id == ctx.guild.owner.id:
        raise CantChoseGuildOwner()
    return True

async def is_owner_block(ctx, user_id: int):
    errors = []
    for check in (is_guild_owner_block,):
        try:
            await check(ctx, user_id)
            return True
        except commands.CheckFailure as e:
            errors.append(str(e))
    raise commands.CheckFailure("\n".join(f"❌ {msg}" for msg in errors))

async def is_owner_developer(ctx):
    errors = []
    for check in (is_developer, is_guild_owner, is_system_owner):
        try:
            if await check(ctx):
                return True
        except commands.CheckFailure as e:
            errors.append(str(e))
    raise commands.CheckFailure("\n".join(f"❌ {msg}" for msg in errors))

async def is_owner_developer_coordinator(ctx):
    errors = []
    for check in (is_coordinator, is_developer, is_guild_owner, is_system_owner):
        try:
            if await check(ctx):
                return True
        except commands.CheckFailure as e:
            errors.append(str(e))
    raise commands.CheckFailure("\n".join(f"❌ {msg}" for msg in errors))

async def is_owner_developer_coordinator_moderator(ctx):
    errors = []
    for check in (is_moderator, is_coordinator, is_developer, is_guild_owner, is_system_owner):
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

async def check_owner_dev_coord_mod(ctx: commands.Context, channel: discord.abc.GuildChannel) -> tuple[bool, bool]:
    is_owner_or_dev = False
    has_any_permission = False
    channel_obj = channel or ctx.channel
    for check in (is_system_owner, is_guild_owner, is_developer):
        try:
            if await check(ctx):
                is_owner_or_dev = True
                has_any_permission = True
                break
        except commands.CheckFailure as e:
            continue
    if not has_any_permission:
        try:
            if await is_coordinator_for_channel(ctx, channel_obj):
                has_any_permission = True
        except commands.CheckFailure as e:
            pass
    if not has_any_permission:
        try:
            if await is_moderator_for_channel(ctx, channel_obj):
                has_any_permission = True
        except commands.CheckFailure as e:
            pass
    return has_any_permission, is_owner_or_dev

async def send_check_failure_embed(ctx: commands.Context, error: commands.CheckFailure, *, title: str = "Permission Check Failed"):
    embed = discord.Embed(
        title=title,
        description=error.args[0] if error.args else "An unknown error occurred.",
        color=discord.Color.red()
    )
    embed.set_footer(text="Contact a bot admin if you think this is a mistake.")
    await ctx.send(embed=embed)
