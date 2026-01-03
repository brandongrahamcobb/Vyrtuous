''' helpers.py The purpose of this program is to provide generic parameters.

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
from os.path import dirname, abspath, expanduser, join

#### DEVELOPMENT
RELEASE_MODE = False
#### DIRECTORIES
DIR_BASE = abspath(join(dirname(dirname(dirname(dirname(__file__))))))
DIR_HOME = expanduser('~')
#### DISCORD
DISCORD_CHARACTER_LIMITS = [2000, 4000]
DISCORD_CHARACTER_LIMIT = 2000
DISCORD_COGS = [
    'vyrtuous.cogs.admin_commands',
    'vyrtuous.cogs.aliases',
    'vyrtuous.cogs.coordinator_commands',
    'vyrtuous.cogs.dev_commands',
    'vyrtuous.cogs.event_listeners',
    'vyrtuous.cogs.everyone_commands',
    'vyrtuous.cogs.guild_owner_commands',
    'vyrtuous.cogs.help_command',
    'vyrtuous.cogs.heartbeat',
    'vyrtuous.cogs.moderator_commands',
    'vyrtuous.cogs.scheduled_tasks',
    'vyrtuous.cogs.system_owner_commands'
]
DISCORD_COGS_CLASSES = [
    'AdminCommands',
    'Aliases',
    'CoordinatorCommands',
    'DevCommands',
    'EventListeners',
    'EveryoneCommands',
    'GuildOwnerCommands',
    'HelpCommand',
    'Heartbeat',
    'ModeratorCommands',
    'ScheduledTasks'
    'SystemOwnerCommands',
]
DISCORD_COMMAND_PREFIX = '!'
#### PATHS
PATH_TOML = join(DIR_HOME, 'git', 'python', 'Vyrtuous', 'pyproject.toml')
PATH_LOG = join(DIR_BASE, 'vyrtuous', '.log', 'discord.log')


def make_mock_snowflake(seed):
    return 10_000_000_000_000_000 + seed

PRIVILEGED_AUTHOR_ID = make_mock_snowflake(1)
PRIVILEGED_AUTHOR_NAME = "BlackBox"

NOT_PRIVILEGED_AUTHOR_ID = make_mock_snowflake(2)
NOT_PRIVILEGED_AUTHOR_NAME = "WhiteBox"

TEXT_CHANNEL_ID = make_mock_snowflake(10)
TEXT_CHANNEL_NAME = "text-channel"

VOICE_CHANNEL_ONE_ID = make_mock_snowflake(11)
VOICE_CHANNEL_ONE_NAME = "Voice Channel 1"

VOICE_CHANNEL_TWO_ID = make_mock_snowflake(12)
VOICE_CHANNEL_TWO_NAME = "Voice Channel 2"

MESSAGE_ID = make_mock_snowflake(100)

ROLE_ID = make_mock_snowflake(200)
ROLE_NAME = "Role"

GUILD_ID = make_mock_snowflake(500)
GUILD_NAME = "Guild"

