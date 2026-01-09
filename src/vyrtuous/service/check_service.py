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
from vyrtuous.enhanced_members.administrator import Administrator
from vyrtuous.enhanced_members.coordinator import Coordinator
from vyrtuous.enhanced_members.developer import Developer
from vyrtuous.enhanced_members.moderator import Moderator
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
    def __init__(self, message='You are not the system owner and cannot do this.'):
        super().__init__(message)

class NotDeveloper(commands.CheckFailure):
    def __init__(self, message='You are not a developer in this server.'):
        super().__init__(message)
                                    
class NotAtHome(commands.CheckFailure):
    def __init__(self, message='You are not in the home Discord!.'):
        super().__init__(message)

def at_home(ctx_interaction_or_message) -> bool:
    bot = DiscordBot.get_instance()
    if ctx_interaction_or_message.guild.id == int(bot.config['discord_testing_guild_snowflake']):
        return True
    raise NotAtHome()

async def is_moderator(ctx_interaction_or_message):
    if isinstance(ctx_interaction_or_message, discord.Interaction):
        member_snowflake = ctx_interaction_or_message.user.id
    elif isinstance(ctx_interaction_or_message, commands.Context):
        member_snowflake = ctx_interaction_or_message.author.id
    elif isinstance(ctx_interaction_or_message, discord.Message):
        member_snowflake = ctx_interaction_or_message.author.id
    else:
        member_snowflake = None
    moderator = await Moderator.fetch_by_guild_and_member(guild_snowflake=ctx_interaction_or_message.guild.id, member_snowflake=member_snowflake)
    if not moderator:
        raise NotModerator()
    return True

async def is_coordinator(ctx_interaction_or_message):
    if isinstance(ctx_interaction_or_message, discord.Interaction):
        member_snowflake = ctx_interaction_or_message.user.id
    elif isinstance(ctx_interaction_or_message, commands.Context):
        member_snowflake = ctx_interaction_or_message.author.id
    elif isinstance(ctx_interaction_or_message, discord.Message):
        member_snowflake = ctx_interaction_or_message.author.id
    else:
        member_snowflake = None
    coordinator = await Coordinator.fetch_by_guild_and_member(guild_snowflake=ctx_interaction_or_message.guild.id, member_snowflake=member_snowflake)
    if not coordinator:
        raise NotCoordinator()
    return True
    
async def is_administrator(ctx_interaction_or_message):
    if isinstance(ctx_interaction_or_message, discord.Interaction):
        member_snowflake = ctx_interaction_or_message.user.id
    elif isinstance(ctx_interaction_or_message, commands.Context):
        member_snowflake = ctx_interaction_or_message.author.id
    elif isinstance(ctx_interaction_or_message, discord.Message):
        member_snowflake = ctx_interaction_or_message.author.id
    else:
        member_snowflake = None
    administrator = await Administrator.fetch_by_guild_and_member(guild_snowflake=ctx_interaction_or_message.guild.id, member_snowflake=member_snowflake)
    if not administrator:
        raise NotAdministrator()
    return True

async def is_developer(ctx_interaction_or_message):
    if isinstance(ctx_interaction_or_message, discord.Interaction):
        member_snowflake = ctx_interaction_or_message.user.id
    elif isinstance(ctx_interaction_or_message, commands.Context):
        member_snowflake = ctx_interaction_or_message.author.id
    elif isinstance(ctx_interaction_or_message, discord.Message):
        member_snowflake = ctx_interaction_or_message.author.id
    else:
        member_snowflake = None
    developer = await Developer.fetch_by_guild_and_member(guild_snowflake=ctx_interaction_or_message.guild.id, member_snowflake=member_snowflake)
    if developer:
        return True
    else:
        raise NotDeveloper()

async def is_guild_owner(ctx_interaction_or_message):
    if isinstance(ctx_interaction_or_message, discord.Interaction):
        member_snowflake = ctx_interaction_or_message.user.id
    elif isinstance(ctx_interaction_or_message, commands.Context):
        member_snowflake = ctx_interaction_or_message.author.id
    elif isinstance(ctx_interaction_or_message, discord.Message):
        member_snowflake = ctx_interaction_or_message.author.id
    else:
        member_snowflake = None
    if ctx_interaction_or_message.guild.owner_id != member_snowflake:
        raise NotGuildOwner()
    return True

async def is_system_owner(ctx_interaction_or_message):
    bot = DiscordBot.get_instance()
    if isinstance(ctx_interaction_or_message, discord.Interaction):
        user = ctx_interaction_or_message.user
    elif isinstance(ctx_interaction_or_message, commands.Context):
        user = ctx_interaction_or_message.author
    elif isinstance(ctx_interaction_or_message, discord.Message):
        user = ctx_interaction_or_message.author
    system_owner_id = int(bot.config['discord_owner_id'])
    if system_owner_id != user.id:
        raise NotSystemOwner()
    return True

def is_system_owner_developer_predicator():
    async def predicate(ctx_interaction_or_message):
        for check in (is_system_owner, is_developer):
            try:
                if await check(ctx_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not a system owner or developer.')
    predicate._permission_level = 'Developer'
    return commands.check(predicate)

def is_system_owner_developer_guild_owner_predicator():
    async def predicate(ctx_interaction_or_message):
        for check in (is_system_owner, is_developer, is_guild_owner):
            try:
                if await check(ctx_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not a system owner, developer or guild owner in this server.')
    predicate._permission_level = 'Guild Owner'
    return commands.check(predicate)
                          
def is_system_owner_predicator():
    async def predicate(ctx_interaction_or_message):
        if await is_system_owner(ctx_interaction_or_message):
            return True
        raise commands.CheckFailure('You are not a system owner.')
    predicate._permission_level = 'System Owner'
    return commands.check(predicate)

def is_system_owner_developer_guild_owner_administrator_predicator():
    async def predicate(ctx_interaction_or_message):
        for check in (is_system_owner, is_developer, is_guild_owner, is_administrator):
            try:
                if await check(ctx_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not a system owner, developer, guild owner or administrator in this server.')
    predicate._permission_level = 'Administrator'
    return commands.check(predicate)
    
def is_system_owner_developer_guild_owner_administrator_coordinator_predicator():
    async def predicate(ctx_interaction_or_message):
        for check in (is_system_owner, is_developer, is_guild_owner, is_administrator, is_coordinator):
            try:
                if await check(ctx_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not a system owner, developer, guild owner, administrator or coordinator in this channel.')
    predicate._permission_level = 'Coordinator'
    return commands.check(predicate)
    
def is_system_owner_developer_guild_owner_administrator_coordinator_moderator_predicator():
    async def predicate(ctx_interaction_or_message):
        for check in (is_system_owner, is_developer, is_guild_owner, is_administrator, is_coordinator, is_moderator):
            try:
                if await check(ctx_interaction_or_message):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure('You are not a system owner, developer, guild owner, administrator, coordinator or moderator in this channel.')
    predicate._permission_level = 'Moderator'
    return commands.check(predicate)

# async def permission_check(ctx_interaction_or_message) -> str:
#     checks = (
#         ('System Owner', is_system_owner),
#         ('Developer', is_developer),
#         ('Guild Owner', is_guild_owner),
#         ('Administrator', is_administrator),
#         ('Coordinator', is_coordinator),
#         ('Moderator', is_moderator),
#     )
#     for role_name, check in checks:
#         try:
#             if await check(ctx_interaction_or_message):
#                 return role_name
#         except commands.CheckFailure:
#             continue
#     return 'Everyone'

async def permission_check(ctx_interaction_or_message, omit=()) -> str:
    checks = (
        ('System Owner', is_system_owner),
        ('Developer', is_developer),
        ('Guild Owner', is_guild_owner),
        ('Administrator', is_administrator),
        ('Coordinator', is_coordinator),
        ('Moderator', is_moderator),
    )
    omit_set = set(omit)
    for role_name, check in checks:
        if role_name in omit_set:
            continue
        try:
            if await check(ctx_interaction_or_message):
                return role_name
        except commands.CheckFailure:
            continue
    return 'Everyone'

async def member_is_guild_owner(guild_snowflake: int, member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild and guild.owner_id == member_snowflake:
        return True
    raise NotGuildOwner

async def member_is_system_owner(member_snowflake: int) -> bool:
    bot = DiscordBot.get_instance()
    if int(bot.config['discord_owner_id']) == member_snowflake:
        return True
    raise NotSystemOwner

async def member_is_developer(guild_snowflake: int, member_snowflake: int) -> bool:
    developer = await Developer.fetch_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
    if not developer:
        raise NotDeveloper
    return True

async def member_is_administrator(guild_snowflake: int, member_snowflake: int) -> bool:
    administrator = await Administrator.fetch_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
    if not administrator:
        raise NotAdministrator
    return True

async def member_is_coordinator(channel_snowflake: int, guild_snowflake: int, member_snowflake: int) -> bool:
    coordinator = await Coordinator.fetch_by_channel_guild_and_member(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
    if not coordinator:
        raise NotCoordinator
    return True

async def member_is_moderator(channel_snowflake: int, guild_snowflake: int, member_snowflake: int) -> bool:
    moderator = await Moderator.fetch_by_channel_guild_and_member(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
    if not moderator:
        raise NotModerator
    return True

async def has_equal_or_higher_role(message_ctx_interaction, channel_snowflake: int, guild_snowflake: int, member_snowflake: int, sender_snowflake: int) -> bool:
    bot = DiscordBot.get_instance() 
    async def get_highest_role(member_sf):
        try:
            if await member_is_system_owner( member_snowflake=member_sf):
                return PERMISSION_TYPES.index('System Owner')
        except NotSystemOwner:
            pass
        try:
            if await member_is_developer(guild_snowflake=guild_snowflake, member_snowflake=member_sf):
                return PERMISSION_TYPES.index('Developer')
        except NotDeveloper:
            pass
        try:
            if await member_is_guild_owner(guild_snowflake=guild_snowflake, member_snowflake=member_sf):
                return PERMISSION_TYPES.index('Guild Owner')
        except NotGuildOwner:
            pass
        try:
            if await member_is_administrator(guild_snowflake=guild_snowflake, member_snowflake=member_sf):
                return PERMISSION_TYPES.index('Administrator')
        except NotAdministrator:
            pass
        if channel_snowflake:
            try:
                if await member_is_coordinator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_sf):
                    return PERMISSION_TYPES.index('Coordinator')
            except NotCoordinator:
                pass
            try:
                if await member_is_moderator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_sf):
                    return PERMISSION_TYPES.index('Moderator')
            except NotModerator:
                pass
        return PERMISSION_TYPES.index('Everyone')
    sender_rank = await get_highest_role(sender_snowflake)
    if isinstance(message_ctx_interaction, (commands.Context, discord.Message)):
        if message_ctx_interaction.author.id == member_snowflake:
            return PERMISSION_TYPES[sender_rank]
    elif isinstance(message_ctx_interaction, discord.Interaction):
        if message_ctx_interaction.user.id == member_snowflake:
            return PERMISSION_TYPES[sender_rank]
    target_rank = await get_highest_role(member_snowflake)
    msg = f'You may not execute this command on this `{PERMISSION_TYPES[target_rank]}` because they have equal or higher role than you in this channel/server.'
    if sender_rank <= target_rank:
        if isinstance(message_ctx_interaction, discord.Interaction):
            raise app_commands.CheckFailure(msg)
        else:
            raise commands.CheckFailure(msg)
    return PERMISSION_TYPES[sender_rank]

async def is_system_owner_developer_guild_owner_administrator_coordinator_via_channel_member(channel_snowflake: int, guild_snowflake: int, member_snowflake: int) -> str:
    checks = (
        ('System Owner', lambda: member_is_system_owner(member_snowflake=member_snowflake)),
        ('Guild Owner', lambda: member_is_guild_owner(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)),
        ('Developer', lambda: member_is_developer(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)),
        ('Administrator', lambda: member_is_administrator(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)),
        ('Coordinator', lambda: member_is_coordinator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake))
    )
    for role_name, check in checks:
        try:
            if await check():
                return role_name
        except (NotSystemOwner, NotGuildOwner, NotDeveloper, NotAdministrator, NotCoordinator):
            continue
    return 'Everyone'

async def is_system_owner_developer_guild_owner_administrator_coordinator_moderator_via_channel_member(channel_snowflake: int, guild_snowflake: int, member_snowflake: int) -> str:
    checks = (
        ('System Owner', lambda: member_is_system_owner(member_snowflake=member_snowflake)),
        ('Developer', lambda: member_is_developer(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)),
        ('Guild Owner', lambda: member_is_guild_owner(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)),
        ('Administrator', lambda: member_is_administrator(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)),
        ('Coordinator', lambda: member_is_coordinator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)),
        ('Moderator', lambda: member_is_moderator(channel_snowflake=channel_snowflake, guild_snowflake=guild_snowflake, member_snowflake=member_snowflake))
    )
    for role_name, check in checks:
        try:
            if await check():
                return role_name
        except (NotSystemOwner, NotGuildOwner, NotDeveloper, NotAdministrator, NotCoordinator, NotModerator):
            continue
    return 'Everyone'

def check_not_self(ctx_interaction_or_message, member_snowflake: int):
    try:
        if member_snowflake == ctx_interaction_or_message.guild.me.id:
            raise commands.CheckFailure('You cannot execute actions on {ctx_interaction_or_message.guild.me.mention}.')
        else:
            return True
    except commands.CheckFailure:
        raise        
    