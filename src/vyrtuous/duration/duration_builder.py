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

from datetime import datetime

from vyrtuous.duration.duration_service import DurationService


class DurationBuilder:
    def __init__(self):
        self.__duration_service = DurationService()
        self.__duration = None

    def parse(self, value):
        self.__duration = self.__duration_service.parse(value)
        return self

    def from_seconds(self, seconds):
        self.__duration = self.__duration_service.from_seconds(seconds)
        return self

    def from_timedelta(self, td, prefix="+"):
        self.__duration = self.__duration_service.from_timedelta(td, prefix)
        return self

    def to_seconds(self):
        return self.__duration_service.to_seconds(self.__duration)

    def to_timedelta(self):
        return self.__duration_service.to_timedelta(self.__duration)

    def to_expires_in(self, base: datetime | None = None):
        return self.__duration_service.to_expires_in(self.__duration, base=base)

    def from_timestamp(self, expires_in: datetime):
        self.__duration = self.__duration_service.from_timestamp(expires_in)
        return self

    def build(self, *, as_str: bool = False):
        if as_str and self.__duration:
            return f"{self.__duration.prefix}{self.__duration.number}{self.__duration.unit}"
        return self.__duration

    def to_unix_ts(self, base: datetime | None = None):
        if self.__duration.number == 0:
            return "permanent"
        return f"<t:{int(self.to_expires_in(base).timestamp())}:R>"
