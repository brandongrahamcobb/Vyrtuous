"""!/bin/python3

time_to_complete.py The purpose of this program is to provide the TimeToComplete utility module.

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


class TimeToComplete:

    def is_around_one_second(self, elapsed: float = 1.0):
        return 0.0 <= elapsed <= 2.0

    def time_elapsed_measurement(self, start: datetime, end: datetime) -> float:
        if start is None or end is None:
            return 0.0
        return (end - start).total_seconds()
