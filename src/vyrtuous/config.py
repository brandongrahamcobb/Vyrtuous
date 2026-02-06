"""!/bin/python3
config.py  The purpose of this program is to provide the primary configuaration script.

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

import os
from typing import Any, Dict

from dotenv import load_dotenv


class Config:

    _config = None

    @classmethod
    def get_config(cls) -> Dict[str, Any]:
        load_dotenv()
        _config = {
            "black_box_api_key": os.environ["BLACK_BOX_API_KEY"],
            "white_box_api_key": os.environ["WHITE_BOX_API_KEY"],
            "vyrtuous_api_key": os.environ["VYRTUOUS_API_KEY"],
            "client_id": os.environ["CLIENT_ID"],
            "client_secret": os.environ["CLIENT_SECRET"],
            "redirect_uri": os.environ["REDIRECT_URI"],
            "discord_command_prefix": os.environ["DISCORD_COMMAND_PREFIX"],
            "discord_owner_id": int(os.environ["DISCORD_OWNER_ID"]),
            "discord_testing_guild_snowflake": int(
                os.environ["DISCORD_TESTING_GUILD_ID"]
            ),
            "logging_level": os.environ["LOGGING_LEVEL"],
            "postgres_database": os.getenv("POSTGRES_DB"),
            "backup_directory": os.getenv("DB_DIRECTORY"),
            "postgres_host": os.getenv("POSTGRES_HOST"),
            "postgres_password": os.getenv("POSTGRES_PASSWORD"),
        }
        _config["release_mode"] = str(os.environ.get("RELEASE_MODE")).lower() in (
            "1",
            "true",
        )
        return _config
