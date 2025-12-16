''' config.py  The purpose of this program is to provide the primary configuaration script.

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
from dotenv import load_dotenv
from typing import Any, Dict
from vyrtuous.inc.helpers import *

import os

class Config:

    _config = None

    @classmethod
    def get_config(cls) -> Dict[str, Any]:
        load_dotenv()
        _config = {
            'discord_api_key': os.environ['DISCORD_API_KEY'],
            'client_id': os.environ['CLIENT_ID'],
            'client_secret': os.environ['CLIENT_SECRET'],
            'redirect_uri': os.environ['REDIRECT_URI'],
            'discord_command_prefix': os.environ['DISCORD_COMMAND_PREFIX'],
            'discord_owner_id': os.environ['DISCORD_OWNER_ID'],
            'discord_testing_guild_id': os.environ['DISCORD_TESTING_GUILD_ID'],
            'logging_level': os.environ['LOGGING_LEVEL'],
            'release_mode': os.environ['RELEASE_MODE']
        }
        return _config
