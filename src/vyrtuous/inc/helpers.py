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

#### DEVELOPMENT
RELEASE_MODE = False
#### DIRECTORIES
DIR_BASE = abspath(join(dirname(dirname(dirname(dirname(__file__))))))
DIR_HOME = expanduser("~")
#### DISCORD
DISCORD_CHARACTER_LIMITS = [2000, 4000]
DISCORD_CHARACTER_LIMIT = 2000
DISCORD_COGS = [
    "vyrtuous.administrator.administrator_app_commands",
    "vyrtuous.administrator.administrator_text_commands",
    "vyrtuous.coordinator.coordinator_app_commands",
    "vyrtuous.coordinator.coordinator_text_commands",
    "vyrtuous.developer.developer_app_commands",
    "vyrtuous.developer.developer_text_commands",
    "vyrtuous.owner.guild_owner_app_commands",
    "vyrtuous.owner.guild_owner_text_commands",
    "vyrtuous.moderator.help_app_command",
    "vyrtuous.moderator.help_text_command",
    "vyrtuous.moderator.moderator_app_commands",
    "vyrtuous.moderator.moderator_text_commands",
    "vyrtuous.sysadmin.sysadmin_app_commands",
    "vyrtuous.sysadmin.sysadmin_text_commands",
    "vyrtuous.utils.channel_event_listeners",
    "vyrtuous.utils.generic_event_listeners",
    "vyrtuous.utils.guild_event_listeners",
    "vyrtuous.utils.heartbeat",
    "vyrtuous.utils.scheduled_tasks",
]
DISCORD_COGS_CLASSES = [
    "AdminAppCommands",
    "AdminTextCommands",
    "ChannelEventListeners",
    "CoordinatorAppCommands",
    "CoordinatorTextCommands",
    "DevAppCommands",
    "DevTextCommands",
    "GenericEventListeners",
    "GuildEventListeners",
    "GuildOwnerAppCommands",
    "GuildOwnerTextCommands",
    "HelpAppCommand",
    "HelpTextCommand",
    "Heartbeat",
    "ModeratorAppCommands",
    "ModeratorTextCommands",
    "ScheduledTasks",
    "SysadminAppCommands",
    "SysadminTextCommands",
]
DISCORD_COMMAND_PREFIX = "!"
#### PATHS
PATH_TOML = join(DIR_HOME, "git", "sandbox", "python", "Vyrtuous", "pyproject.toml")
PATH_LOG = join(DIR_BASE, "vyrtuous", ".log", "discord.log")
