"""!/bin/python3
helpers.py The purpose of this program is to provide generic parameters.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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
"""

from os.path import abspath, dirname, expanduser, join
from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot

#### DEVELOPMENT
RELEASE_MODE = False
#### DIRECTORIES
DIR_BASE = abspath(join(dirname(dirname(dirname(dirname(__file__))))))
DIR_HOME = expanduser("~")
#### DISCORD
DISCORD_CHARACTER_LIMITS = [2000, 4000]
DISCORD_CHARACTER_LIMIT = 2000
DISCORD_COGS = [
    "vyrtuous.app_commands.administrator_app_commands",
    "vyrtuous.app_commands.coordinator_app_commands",
    "vyrtuous.app_commands.developer_app_commands",
    "vyrtuous.app_commands.guild_owner_app_commands",
    "vyrtuous.app_commands.help_app_command",
    "vyrtuous.app_commands.hidden_administrator_app_commands",
    "vyrtuous.app_commands.hidden_coordinator_app_commands",
    "vyrtuous.app_commands.hidden_developer_app_commands",
    "vyrtuous.app_commands.hidden_guild_owner_app_commands",
    "vyrtuous.app_commands.hidden_moderator_app_commands",
    "vyrtuous.app_commands.hidden_sysadmin_app_commands",
    "vyrtuous.app_commands.moderator_app_commands",
    "vyrtuous.app_commands.sysadmin_app_commands",
    "vyrtuous.listeners.channel_event_listeners",
    "vyrtuous.listeners.generic_event_listeners",
    "vyrtuous.listeners.guild_event_listeners",
    "vyrtuous.listeners.scheduled_tasks",
    "vyrtuous.listeners.startup_listener",
    "vyrtuous.system.heartbeat",
    "vyrtuous.text_commands.administrator_text_commands",
    "vyrtuous.text_commands.coordinator_text_commands",
    "vyrtuous.text_commands.developer_text_commands",
    "vyrtuous.text_commands.guild_owner_text_commands",
    "vyrtuous.text_commands.help_text_command",
    "vyrtuous.text_commands.hidden_administrator_text_commands",
    "vyrtuous.text_commands.hidden_coordinator_text_commands",
    "vyrtuous.text_commands.hidden_developer_text_commands",
    "vyrtuous.text_commands.hidden_guild_owner_text_commands",
    "vyrtuous.text_commands.hidden_moderator_text_commands",
    "vyrtuous.text_commands.hidden_sysadmin_text_commands",
    "vyrtuous.text_commands.moderator_text_commands",
    "vyrtuous.text_commands.sysadmin_text_commands",
]
DISCORD_COGS_CLASSES = [
    "AdministratorAppCommands",
    "AdministratorTextCommands",
    "ChannelEventListeners",
    "CoordinatorAppCommands",
    "CoordinatorTextCommands",
    "DeveloperAppCommands",
    "DeveloperTextCommands",
    "GenericEventListeners",
    "GuildEventListeners",
    "GuildOwnerAppCommands",
    "GuildOwnerTextCommands",
    "HelpAppCommand",
    "HelpTextCommand",
    "Heartbeat",
    "HiddenAdministratorAppCommands",
    "HiddenAdministratorTextCommands",
    "HiddenCoordinatorAppCommands",
    "HiddenCoordinatorTextCommands",
    "HiddenDeveloperAppCommands",
    "HiddenDeveloperTextCommands",
    "HiddenGuildOwnerAppCommands",
    "HiddenGuildOwnerTextCommands",
    "HiddenModeratorAppCommands",
    "HiddenModeratorTextCommands",
    "HiddenSysadminAppCommands",
    "HiddenSysadminTextCommands",
    "ModeratorAppCommands",
    "ModeratorTextCommands",
    "ScheduledTasks",
    "Startup",
    "SysadminAppCommands",
    "SysadminTextCommands",
]
DISCORD_COMMAND_PREFIX = "!"
#### PATHS
PATH_LOG = join(DIR_BASE, "vyrtuous", ".log", "discord.log")


def at_home(
    source: Union[commands.Context, discord.Interaction, discord.Message],
) -> bool:
    bot: DiscordBot = DiscordBot.get_instance()
    if source.guild is None:
        return False
    if source.guild.id == int(bot.config["discord_testing_guild_snowflake"]):
        return True
    return False


def resolve_author(source) -> discord.User | discord.Member:
    if isinstance(source, discord.Interaction):
        member = source.user
    elif isinstance(source, (commands.Context, discord.Message)):
        member = source.author
    else:
        raise commands.MemberNotFound("Source")
    return member
