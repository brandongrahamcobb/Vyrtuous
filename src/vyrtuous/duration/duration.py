"""!/bin/python3
duration.py The purpose of this program is to provide the Duration properties class.

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

from datetime import datetime, timezone


class Duration:
    def __init__(
        self,
        number: int,
        prefix: str,
        sign: int,
        unit: str,
        *,
        base: datetime = datetime.now(timezone.utc),
    ):
        self.base = base
        self.number = number
        self.prefix = prefix
        self.sign = sign
        self.unit = unit

    # def __str__(self):
    #     if self._number == 0:
    #         return "permanent"
    #     return f"<t:{int(self.expires_in.timestamp())}:R>"
