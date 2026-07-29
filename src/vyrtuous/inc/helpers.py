"""!/bin/python3
helpers.py The purpose of this program is to provide generic parameters and functions.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from os.path import abspath, dirname, expanduser, join

import discord
from discord.ext import commands

from vyrtuous.utils.errors.error import MemberNotFound

#### DEVELOPMENT
RELEASE_MODE = False
#### DIRECTORIES
DIR_BASE = abspath(join(dirname(dirname(dirname(dirname(__file__))))))
DIR_HOME = expanduser("~")
#### DISCORD
DISCORD_CHARACTER_LIMITS = [2000, 4000]
DISCORD_CHARACTER_LIMIT = 2000
DISCORD_COGS = [
    "vyrtuous.app_commands.channel_management_app_commands",
    "vyrtuous.app_commands.development_app_commands",
    # "vyrtuous.app_commands.help_app_command",
    "vyrtuous.app_commands.info_app_commands",
    "vyrtuous.app_commands.moderation_app_commands",
    "vyrtuous.app_commands.user_management_app_commands",
    "vyrtuous.app_commands.utility_app_commands",
    "vyrtuous.listeners.channel_event_listeners",
    "vyrtuous.listeners.generic_event_listeners",
    "vyrtuous.listeners.guild_event_listeners",
    "vyrtuous.listeners.scheduled_tasks",
    "vyrtuous.listeners.startup_listener",
    "vyrtuous.system.heartbeat",
    "vyrtuous.text_commands.alias_management_text_commands",
    "vyrtuous.text_commands.channel_management_text_commands",
    "vyrtuous.text_commands.development_text_commands",
    # "vyrtuous.text_commands.help_text_command",
    "vyrtuous.text_commands.info_text_commands",
    "vyrtuous.text_commands.moderation_text_commands",
    "vyrtuous.text_commands.user_management_text_commands",
    "vyrtuous.text_commands.utility_text_commands",
]
DISCORD_COGS_CLASSES = [
    "ChannelManagementAppCommands",
    "ChannelManagementTextCommands",
    "ChannelEventListeners",
    "DevelopmentAppCommands",
    "DevelopmentTextCommands",
    "GenericEventListeners",
    "GuildEventListeners",
    "Heartbeat",
    "HelpAppCommand",
    "HelpTextCommand",
    "InfoAppCommands",
    "InfoTextCommands",
    "ModerationAppCommands",
    "ModerationTextCommands",
    "UserManagementAppCommands",
    "UserManagementTextCommands",
    "UtilityAppCommands",
    "UtilityTextCommands",
    "ScheduledTasks",
    "Startup",
]
DISCORD_COMMAND_PREFIX = "!"
#### PATHS
PATH_LOG = join(DIR_BASE, ".log", "discord.log")
PATH_GROUPS = join(DIR_BASE, "groups.yml")


def resolve_author(source) -> discord.User | discord.Member:
    if isinstance(source, discord.Interaction):
        member = source.user
    elif isinstance(source, (commands.Context, discord.Message)):
        member = source.author
    else:
        raise MemberNotFound("Source")
    return member
