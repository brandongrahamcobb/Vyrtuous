''' helpers.py  The purpose of this program is to provide generic parameters.
    Copyright (C) 2024 github.com/brandongrahamcobb

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

import discord
from vyrtuous.utils.inc.setup_logging import logger


def parse_comma_number(s):
    return int(s.replace(",", ""))

#### DIRECTORIES
DIR_BASE = abspath(join(dirname(dirname(dirname(dirname(__file__))))))
DIR_HOME = expanduser('~')
DIR_TEMP = join(DIR_BASE, 'vyrtuous', 'temp')

#### PATHS
### Path: Root (Py_vyrtuous)
PATH_TOML = join(DIR_HOME, 'git', 'python', 'Vyrtuous', 'pyproject.toml')
## Path: Source (src)
# Path: Py_Py_vyrtuous (vyrtuous)
PATH_BASHRC = join(DIR_HOME, '.bashrc_vyrtuous')
PATH_CONFIG = join(DIR_BASE, 'vyrtuous', 'config.py')
PATH_CONFIG_YAML = join(DIR_HOME, '.config', 'vyrtuous', 'config.yaml')
PATH_LOG = join(DIR_BASE, 'vyrtuous', '.log', 'discord.log')
PATH_MAIN = join(DIR_BASE, 'vyrtuous', 'main.py')
PATH_USERS = join(DIR_BASE, 'vyrtuous', '.users', 'users.yaml')
# Path: Bots (bots)
PATH_DISCORD_BOT = join(DIR_BASE, 'vyrtuous', 'bots', 'discord_bot.py')
# Path: Cogs (cogs)
PATH_EVENT_LISTENERS = join(DIR_BASE, 'vyrtuous', 'cogs', 'event_listeners.py')
PATH_OWNER_COMMANDS = join(DIR_BASE, 'vyrtuous', 'cogs', 'commands_extra.py')
PATH_PUBLIC_COMMANDS = join(DIR_BASE, 'vyrtuous', 'cogs', 'commands.py')
PATH_SCHEDULED_TASKS = join (DIR_BASE, 'vyrtuous', 'cogs', 'scheduled_tasks.py')
# Path: Handlers (handlers)
PATH_MESSAGE = join(DIR_BASE, 'vyrtuous', 'utils', 'handlers', 'message_service..py')
PATH_PREDICATOR = join(DIR_BASE, 'vyrtuous', 'utils', 'handlers', 'predicator.py')
PATH_ROLE_MANAGER = join(DIR_BASE, 'vyrtuous', 'utils', 'handlers', 'mute_service.py')
# Path: Include (inc)
PATH_HELPERS = join(DIR_BASE, 'vyrtuous', 'utils', 'inc', 'helpers.py')
PATH_INCREMENT_VERSION = join(DIR_BASE, 'vyrtuous', 'utils', 'inc', 'increment_version.py')
PATH_LOAD_CONTENTS = join(DIR_BASE, 'vyrtuous', 'utils', 'inc', 'load_contents.py')
PATH_LOAD_YAML = join(DIR_BASE, 'vyrtuous', 'utils', 'inc', 'load_yaml.py')
PATH_PROMPT_FOR_VALUES = join(DIR_BASE, 'vyrtuous', 'utils', 'inc', 'prompt_for_values.py')
PATH_SETUP_LOGGING = join(DIR_BASE, 'vyrtuous', 'utils', 'inc', 'setup_logging.py')
#### CONTENTS
#### DISCORD
DISCORD_CHARACTER_LIMITS = [2000, 4000]
DISCORD_CHARACTER_LIMIT = 2000
DISCORD_COGS = [
    'vyrtuous.cogs.commands',
    'vyrtuous.cogs.commands_extra',
    'vyrtuous.cogs.event_listeners',
    'vyrtuous.cogs.scheduled_tasks',
]
DISCORD_COMMAND_PREFIX = '!'
DISCORD_DEVELOPER_CHANNEL = 1391574892826329261
DISCORD_INTENTS = discord.Intents.all()
DISCORD_OWNER_ID = 154749533429956608
DISCORD_RELEASE_MODE = False
DISCORD_TESTING_GUILD_ID = 1300517536001036348

LOGGING_LEVEL = 'INFO'

#### SCRIPTURE
SCRIPTURE_HEADERS = {
    'User-Agent': 'brandongrahamcobb@icloud.com',
    'api-key': '2eb327f99245cd3d68da55370656d6e2'
}

USER_AGENT = 'https://github.com/brandongrahamcobb/Py_vyrtuous.git'
VERSION = '1.0.0'
