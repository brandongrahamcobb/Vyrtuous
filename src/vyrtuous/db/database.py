"""database.py The purpose of this program is to provide the database utility module.

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

from datetime import datetime
import os
import subprocess

import asyncpg


class Database:

    def __init__(self, config, *, directory=None):
        self.database: str = config.get("postgres_database")
        self.directory = directory if directory else config.get("backup_directory")
        self.host: str = config.get("postgres_host")
        self.password: str = config.get("postgres_password")
        self.timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.file_name: str = os.path.join(
            self.directory, f"backup_{self.timestamp}.sql"
        )
        self.user: str = str(os.getenv("POSTGRES_USER"))

    def create_backup_directory(self):
        os.makedirs(self.directory, exist_ok=True)
        return

    async def database_init(self):
        return await asyncpg.create_pool(
            host=self.host,
            database=self.database,
            user=self.user,
            password=self.password,
            command_timeout=30,
        )

    def execute_backup(self):
        dump_command = [
            "pg_dump",
            "-U",
            self.user,
            "-h",
            self.host,
            "-d",
            self.database,
            "-F",
            "p",
            "-f",
            self.file_name,
        ]
        env = os.environ.copy()
        env["PGPASSWORD"] = self.password
        result = subprocess.run(
            dump_command,
            capture_output=True,
            text=True,
            env=env,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Backup failed: {result.stderr}")
        return result
