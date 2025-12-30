''' check_service.py The purpose of this program is to provide the check_service module.

    Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
from discord import app_commands
from discord.ext import commands
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.coordinator import Coordinator
from vyrtuous.utils.developer import Developer
from vyrtuous.utils.moderator import Moderator
from vyrtuous.utils.permission import PERMISSION_TYPES
from vyrtuous.utils.setup_logging import logger

import discord
      
class NotCoordinator(commands.CheckFailure):
    def __init__(self, message='You are not a coordinator in the requested channel.'):
        super().__init__(message)

class NotAdministrator(commands.CheckFailure):
    def __init__(self, message='You are not an administrator in the requested channel.'):
        super().__init__(message)
        
class NotModerator(commands.CheckFailure):
    def __init__(self, message='You are not a moderator in the requested channel.'):
        super().__init__(message)

class NotGuildOwner(commands.CheckFailure):
    def __init__(self, message='You are not the guild owner and cannot do this.'):
        super().__init__(message)

class NotSystemOwner(commands.CheckFailure):
    def __init__(self, message='You are not the bot owner.'):
        super().__init__(message)

class NotDeveloper(commands.CheckFailure):
    def __init__(self, message='You are not a developer in this guild.'):
        super().__init__(message)
                                    
class NotAtHome(commands.CheckFailure):
    def __init__(self, message='You are not in the home Discord!.'):
        super().__init__(message)

async def at_home(ctx_or_interaction_or_message) -> bool:
    bot = DiscordBot.get_instance()
    if ctx_or_interaction_or_message.guild.id == int(bot.config['discord_testing_guild_id']):
        return True
    raise NotAtHome()
    
async def send_check_failure_embed(ctx_or_interaction_or_message, error: commands.CheckFailure, *, title: str = 'Permission Check Failed'):
    embed = discord.Embed(
        title=title,
        description=error.args[0] if error.args else 'An unknown error occurred.',
        color=discord.Color.red()
    )
    embed.set_footer(text='Contact a bot admin if you think this is a mistake.')
    if isinstance(ctx_or_interaction_or_message, discord.Interaction):
        return await ctx_or_interaction_or_message.response.send_message(embed=embed)
    elif isinstance(ctx_or_interaction_or_message, commands.Context) or isinstance(ctx_or_interaction_or_message, discord.Message):
        return await ctx_or_interaction_or_message.reply(embed=embed)
    else:
        return None

async def is_moderator(ctx_or_interaction_or_message):
    bot = DiscordBot.get_instance()
    if isinstance(ctx_or_interaction_or_message, discord.Interaction):
        user_id = ctx_or_interaction_or_message.user.id
    elif isinstance(ctx_or_interaction_or_message, commands.Context):
        user_id = ctx_or_interaction_or_message.author.id
    elif isinstance(ctx_or_interaction_or_message, discord.Message):
        user_id = ctx_or_interaction_or_message.author.id
    else:
        user_id = None
    moderator_channel_ids = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=ctx_or_interaction_or_message.guild.id, member_snowflake=user_id)
    if not moderator_channel_ids:
        raise NotModerator()
    return True

async def is_coordinator(ctx_or_interaction_or_message):
    bot = DiscordBot.get_instance()
    if isinstance(ctx_or_interaction_or_message, discord.Interaction):
        user_id = ctx_or_interaction_or_message.user.id
    elif isinstance(ctx_or_interaction_or_message, commands.Context):
        user_id = ctx_or_interaction_or_message.author.id
    elif isinstance(ctx_or_interaction_or_message, discord.Message):
        user_id = ctx_or_interaction_or_message.author.id
    else:
        user_id = None
    channel_snowflakes = await Coordinator.fetch_channels_by_guild_and_member(guild_snowflake=ctx_or_interaction_or_message.guild.id, member_snowflake=user_id)
    if not channel_snowflakes:
        raise NotCoordinator()
    return True
    
async def is_administrator(ctx_or_interaction_or_message):
    if isinstance(ctx_or_interaction_or_message, discord.Interaction):
        user_id = ctx_or_interaction_or_message.user.id
    elif isinstance(ctx_or_interaction_or_message, commands.Context):
        user_id = ctx_or_interaction_or_message.author.id
    elif isinstance(ctx_or_interaction_or_message, discord.Message):
        user_id = ctx_or_interaction_or_message.author.id
    else:
        user_id = None
    administrator = await Administrator.fetch_member(user_id)
    if not administrator:
        raise NotAdministrator()
    return True

async def is_developer(ctx_or_interaction_or_message):
    bot = DiscordBot.get_instance()
    if isinstance(ctx_or_interaction_or_message, discord.Interaction):
        user_id = ctx_or_interaction_or_message.user.id
    elif isinstance(ctx_or_interaction_or_message, commands.Context):
        user_id = ctx_or_interaction_or_message.author.id
    elif isinstance(ctx_or_interaction_or_message, discord.Message):
        user_id = ctx_or_interaction_or_message.author.id
    else:
        user_id = None
    developers = await Developer.fetch_members_by_guild(guild_snowflake=ctx_or_interaction_or_message.guild.id)
    for developer in developers:
        if developer.member_snowflake == user_id:
            found = True
    if not found:
        raise NotDeveloper()
    return True

async def is_guild_owner(ctx_or_interaction_or_message):
    if isinstance(ctx_or_interaction_or_message, discord.Interaction):
        user_id = ctx_or_interaction_or_message.user.id
    elif isinstance(ctx_or_interaction_or_message, commands.Context):
        user_id = ctx_or_interaction_or_message.author.id
    elif isinstance(ctx_or_interaction_or_message, discord.Message):
        user_id = ctx_or_interaction_or_message.author.id
    else:
        user_id = None
    if ctx_or_interaction_or_message.guild.owner_id != user_id:
        raise NotGuildOwner()
    return True

async def is_system_owner(ctx_or_interaction_or_message):
    bot = DiscordBot.get_instance()
    if isinstance(ctx_or_interaction_or_message, discord.Interaction):
        user = ctx_or_interaction_or_message.user
    elif isinstance(ctx_or_interaction_or_message, commands.Context):
        user = ctx_or_interaction_or_message.author
    elif isinstance(ctx_or_interaction_or_message, discord.Message):
        user = ctx_or_interaction_or_message.author
    system_owner_id = int(bot.config['discord_owner_id'])
    if system_owner_id != user.id:
        raise NotSystemOwner()
    return True
    
async def is_owner(ctx_or_interaction_or_message):
    try:
        if await is_system_owner(ctx_or_interaction_or_message):
            return True
    except commands.CheckFailure:
        pass
    try:
        if await is_guild_owner(ctx_or_interaction_or_message):
            return True
    except commands.CheckFailure:
        pass
    return False

def is_owner_predicator():
    async def predicate(ctx_or_interaction_or_message):
        if await is_owner(ctx_or_interaction_or_message):
            return True
        raise commands.CheckFailure('You are not an owner in this guild')
    predicate._permission_level = 'Owner'
    return commands.check(predicate)

def is_owner_developer_predicator():
    async def predicate(ctx_or_interaction_or_message):
        for check in (is_owner, is_developer):
            try:
                if await check(ctx_or_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not an owner or developer in this guild.')
    predicate._permission_level = 'Developer'
    return commands.check(predicate)

def is_owner_developer_administrator_predicator():
    async def predicate(ctx_or_interaction_or_message):
        for check in (is_owner, is_developer, is_administrator):
            try:
                if await check(ctx_or_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not an owner, developer or administrator in this guild.')
    predicate._permission_level = 'Administrator'
    return commands.check(predicate)
    
def is_owner_developer_administrator_coordinator_predicator():
    async def predicate(ctx_or_interaction_or_message):
        for check in (is_owner, is_developer, is_administrator, is_coordinator):
            try:
                if await check(ctx_or_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not an owner, developer, administrator in this guild or a coordinator in this channel.')
    predicate._permission_level = 'Coordinator'
    return commands.check(predicate)
    
def is_owner_developer_administrator_coordinator_moderator_predicator():
    async def predicate(ctx_or_interaction_or_message):
        for check in (is_owner, is_developer, is_administrator, is_coordinator, is_moderator):
            try:
                if await check(ctx_or_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not an owner, developer, administrator in this guild or a coordinator/moderator in this channel.')
    predicate._permission_level = 'Moderator'
    return commands.check(predicate)

async def is_owner_developer_administrator_coordinator_moderator(ctx_or_interaction_or_message) -> str:
    checks = (
        ("Owner", is_system_owner),
        ("Owner", is_guild_owner),
        ("Developer", is_developer),
        ("Administrator", is_administrator),
        ("Coordinator", is_coordinator),
        ("Moderator", is_moderator),
    )
    for role_name, check in checks:
        try:
            if await check(ctx_or_interaction_or_message):
                return role_name
        except commands.CheckFailure:
            continue
    return "Everyone"

async def member_is_owner(guild_snowflake: int, member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild.owner_id == member_snowflake:
        return True
    if int(bot.config['discord_owner_id']) == member_snowflake:
        return True
    return False

async def member_is_developer(guild_snowflake: int, member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        row = await conn.fetchrow('SELECT developer_guild_ids FROM users WHERE discord_snowflake=$1', member_snowflake)
    if not row:
        return False
    return guild_snowflake in (row['developer_guild_ids'] or [])

async def member_is_administrator(guild_snowflake: int, member_snowflake: int) -> bool:
    administrator = await Administrator.fetch_member(member_snowflake)
    if not administrator:
        return False
    return guild_snowflake == administrator.guild_snowflake

async def member_is_coordinator(channel_snowflake: int, guild_snowflake: int, member_snowflake: int) -> bool:
    channel_snowflakes = await Coordinator.fetch_channels_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
    if not channel_snowflakes:
        return False
    return channel_snowflake in channel_snowflakes

async def member_is_moderator(channel_snowflake: int, guild_snowflake: int, member_snowflake: int) -> bool:
    channel_snowflakes = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
    if not channel_snowflakes:
        return False
    return channel_snowflake in channel_snowflakes

async def has_equal_or_higher_role(message_ctx_or_interaction, channel_snowflake: int, guild_snowflake: int, member_snowflake: int, sender_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    async def get_highest_role(member_sf):
        if await member_is_owner(guild_snowflake=guild_snowflake, member_snowflake=member_sf):
            return PERMISSION_TYPES.index("Owner")
        elif await member_is_administrator(guild_snowflake=guild_snowflake, member_snowflake=member_sf):
            return PERMISSION_TYPES.index("Developer")
        elif await member_is_developer(guild_snowflake=guild_snowflake, member_snowflake=member_sf):
            return PERMISSION_TYPES.index("Administrator")
        if channel_snowflake:
            if await member_is_coordinator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_sf):
                return PERMISSION_TYPES.index("Coordinator")
            elif await member_is_moderator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_sf):
                return PERMISSION_TYPES.index("Moderator")
        return PERMISSION_TYPES.index("Everyone")
    sender_rank = await get_highest_role(sender_snowflake)
    target_rank = await get_highest_role(member_snowflake)
    msg = f"You may not execute this command on this `{PERMISSION_TYPES[target_rank]}` because they have equal or higher role than you in this channel/server."
    if sender_rank <= target_rank:
        if isinstance(message_ctx_or_interaction, discord.Interaction):
            raise app_commands.CheckFailure(msg)
        else:
            raise commands.CheckFailure(msg)
    return PERMISSION_TYPES[sender_rank]

async def is_owner_developer_administrator_coordinator_via_channel_member(channel_snowflake: int, guild_snowflake: int, member_snowflake: int) -> str:
    checks = (("Administrator", lambda: member_is_administrator(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)), ("Coordinator", lambda: member_is_coordinator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)), ("Developer", lambda: member_is_developer(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)), ("Owner", lambda: member_is_owner(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)))
    for role_name, check in checks:
        try:
            if await check():
                return role_name
        except commands.CheckFailure:
            continue
    return "Everyone"

async def is_owner_developer_administrator_coordinator_moderator_via_channel_member(channel_snowflake: int, guild_snowflake: int, member_snowflake: int) -> str:
    checks = (("Owner", lambda: member_is_owner(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)), ("Developer", lambda: member_is_developer(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)), ("Administrator", lambda: member_is_administrator(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)), ("Coordinator", lambda: member_is_coordinator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)), ("Moderator", lambda: member_is_moderator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)))
    for role_name, check in checks:
        try:
            if await check():
                return role_name
        except commands.CheckFailure:
            continue
    return "Everyone"

def check_not_self(ctx_interaction_or_message, member_snowflake: int):
    try:
        if member_snowflake == ctx_interaction_or_message.guild.me.id:
            raise commands.CheckFailure("You cannot execute actions on {ctx_interaction_or_message.guild.me.mention}.")
        else:
            return True
    except commands.CheckFailure:
        raise        
    