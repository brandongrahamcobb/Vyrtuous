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

from datetime import datetime, timedelta, timezone

from discord.ext import commands

from vyrtuous.duration.duration import Duration


class DurationService:
    DAYS_PER_WEEK = 7
    DAYS_PER_YEAR = 365
    YEAR_UNITS = {"y", "year", "years"}
    WEEK_UNITS = {"w", "week", "weeks"}
    DAY_UNITS = {"d", "day", "days"}
    HOUR_UNITS = {"h", "hr", "hrs", "hour", "hours"}
    MINUTE_UNITS = {"m", "min", "mins", "minute", "minutes"}
    SECOND_UNITS = {"s", "sec", "secs", "second", "seconds"}
    PREFIXES = {"+", "-", "="}

    UNIT_MAP = {
        **dict.fromkeys(YEAR_UNITS, "y"),
        **dict.fromkeys(WEEK_UNITS, "w"),
        **dict.fromkeys(DAY_UNITS, "d"),
        **dict.fromkeys(HOUR_UNITS, "h"),
        **dict.fromkeys(MINUTE_UNITS, "m"),
        **dict.fromkeys(SECOND_UNITS, "s"),
    }

    def from_timedelta(self, td: timedelta, prefix: str = "+") -> Duration:
        total_seconds = int(td.total_seconds())
        year_seconds = self.DAYS_PER_YEAR * 86400
        week_seconds = self.DAYS_PER_WEEK * 86400
        if total_seconds % year_seconds == 0:
            number, suffix = total_seconds // year_seconds, "y"
        elif total_seconds % week_seconds == 0:
            number, suffix = total_seconds // week_seconds, "w"
        elif total_seconds % 86400 == 0:
            number, suffix = total_seconds // 86400, "d"
        elif total_seconds % 3600 == 0:
            number, suffix = total_seconds // 3600, "h"
        elif total_seconds % 60 == 0:
            number, suffix = total_seconds // 60, "m"
        else:
            number, suffix = total_seconds, "s"
        return self.parse(f"{prefix}{number}{suffix}")

    def from_expires_in(self, expires_in: datetime) -> Duration:
        if expires_in is None:
            return Duration(number=0, prefix="", sign="1", unit="")
        now = datetime.now(timezone.utc)
        remaining = expires_in - now
        total_seconds = int(remaining.total_seconds())
        if total_seconds < 0:
            total_seconds = 0
        year_seconds = self.DAYS_PER_YEAR * 86400
        week_seconds = self.DAYS_PER_WEEK * 86400
        if total_seconds % year_seconds == 0:
            number, unit = total_seconds // year_seconds, "y"
        elif total_seconds % week_seconds == 0:
            number, unit = total_seconds // week_seconds, "w"
        elif total_seconds % 86400 == 0:
            number, unit = total_seconds // 86400, "d"
        elif total_seconds % 3600 == 0:
            number, unit = total_seconds // 3600, "h"
        elif total_seconds % 60 == 0:
            number, unit = total_seconds // 60, "m"
        else:
            number, unit = total_seconds, "s"
        return self.parse(f"+{number}{unit}")

    def from_expires_in_to_str(self, expires_in: datetime) -> str:
        if expires_in is None:
            return 0
        now = datetime.now(timezone.utc)
        remaining = expires_in - now
        total_seconds = int(remaining.total_seconds())
        if total_seconds < 0:
            total_seconds = 0
        if total_seconds % 86400 == 0:
            number, unit = total_seconds // 86400, "d"
        elif total_seconds % 3600 == 0:
            number, unit = total_seconds // 3600, "h"
        elif total_seconds % 60 == 0:
            number, unit = total_seconds // 60, "m"
        else:
            number, unit = total_seconds, "s"
        return f"+{number}{unit}"

    def from_seconds(self, seconds: int):
        return self.parse(f"{seconds}s")

    def to_timedelta(self, duration: Duration) -> timedelta:
        match duration.unit:
            case "y":
                return timedelta(
                    days=duration.number * self.DAYS_PER_YEAR * duration.sign
                )
            case "w":
                return timedelta(
                    days=duration.number * self.DAYS_PER_WEEK * duration.sign
                )
            case "d":
                return timedelta(days=duration.number * duration.sign)
            case "h":
                return timedelta(hours=duration.number * duration.sign)
            case "m":
                return timedelta(minutes=duration.number * duration.sign)
            case "s":
                return timedelta(seconds=duration.number * duration.sign)
            case _:
                raise ValueError(f"Unsupported unit: {duration.unit}")

    def to_expires_in(self, duration, *, base: datetime = None) -> datetime:
        base = base or datetime.now(timezone.utc)
        return base + self.to_timedelta(duration=duration)

    def parse(self, value):
        s = value.lower().strip()
        if s == "0":
            number = 0
            unit = "h"
            prefix = ""
            sign = 1
            return Duration(number=number, prefix=prefix, sign=sign, unit=unit)
        if s[0] in "+-":
            sign = 1 if s[0] == "+" else -1
            s = s[1:]
        else:
            sign = 1
        if s and s[0] in self.PREFIXES:
            prefix = s[0]
            s = s[1:]
        else:
            prefix = ""
        num_str = ""
        for char in s:
            if char.isdigit():
                num_str += char
            else:
                break
        if not num_str:
            raise commands.BadArgument(f"No numeric duration found in '{value}'")
        number = int(num_str)
        s = s[len(num_str) :].strip()
        if not s:
            unit = "h"
        else:
            unit = self.UNIT_MAP.get(s, None)
            if not unit:
                for known in self.UNIT_MAP.keys():
                    if s.startswith(known):
                        unit = self.UNIT_MAP[known]
                        break
            if not unit:
                raise commands.BadArgument(f"Invalid duration unit in '{value}'")
        return Duration(number=number, unit=unit, prefix=prefix, sign=sign)
