"""!/bin/python3
mock_database.py The purpose of this program is to support integration testing for Vyrtuous.

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

database = os.getenv("POSTGRES_DB")
host = os.getenv("POSTGRES_HOST")
password = os.getenv("POSTGRES_PASSWORD")
user = os.getenv("POSTGRES_USER")

dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"
