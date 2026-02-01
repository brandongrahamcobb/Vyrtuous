"""helpers.py The purpose of this program is to provide generic parameters.

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

from os.path import dirname, abspath, expanduser, join

#### DEVELOPMENT
RELEASE_MODE = False
#### DIRECTORIES
DIR_BASE = abspath(join(dirname(dirname(dirname(dirname(__file__))))))
DIR_HOME = expanduser("~")
#### DISCORD
DISCORD_CHARACTER_LIMITS = [2000, 4000]
DISCORD_CHARACTER_LIMIT = 2000
DISCORD_COGS = [
    "vyrtuous.cogs.app_commands.admin_app_commands",
    "vyrtuous.cogs.app_commands.coordinator_app_commands",
    "vyrtuous.cogs.app_commands.dev_app_commands",
    "vyrtuous.cogs.app_commands.guild_owner_app_commands",
    "vyrtuous.cogs.app_commands.moderator_app_commands",
    "vyrtuous.cogs.app_commands.sysadmin_app_commands",
    "vyrtuous.cogs.event_listeners.channel_event_listeners",
    "vyrtuous.cogs.event_listeners.generic_event_listeners",
    "vyrtuous.cogs.event_listeners.guild_event_listeners",
    "vyrtuous.cogs.help_command",
    "vyrtuous.cogs.heartbeat",
    "vyrtuous.cogs.text_commands.admin_text_commands",
    "vyrtuous.cogs.text_commands.coordinator_text_commands",
    "vyrtuous.cogs.text_commands.dev_text_commands",
    "vyrtuous.cogs.text_commands.guild_owner_text_commands",
    "vyrtuous.cogs.text_commands.moderator_text_commands",
    "vyrtuous.cogs.text_commands.sysadmin_text_commands",
    "vyrtuous.cogs.scheduled_tasks",
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
    "HelpCommand",
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

PERMISSION_TYPES = [
    "Everyone",
    "Moderator",
    "Coordinator",
    "Administrator",
    "Guild Owner",
    "Developer",
    "Sysadmin",
]
TARGET_PERMISSIONS = (
    "add_reactions",
    "manage_messages",
    "move_members",
    "mute_members",
    "send_messages",
    "view_channel",
)
CHUNK_SIZE = 7
