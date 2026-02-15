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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.duration.duration import DurationObject


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

    @classmethod
    def from_timedelta(cls, td: timedelta, prefix: str = "+") -> "DurationObject":
        total_seconds = int(td.total_seconds())
        year_seconds = cls.DAYS_PER_YEAR * 86400
        week_seconds = cls.DAYS_PER_WEEK * 86400
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
        return DurationObject(f"{prefix}{number}{suffix}")

    @classmethod
    def from_expires_in(cls, expires_in: datetime) -> "DurationObject":
        if expires_in is None:
            return DurationObject("0")
        now = datetime.now(timezone.utc)
        remaining = expires_in - now
        total_seconds = int(remaining.total_seconds())
        if total_seconds < 0:
            total_seconds = 0
        year_seconds = cls.DAYS_PER_YEAR * 86400
        week_seconds = cls.DAYS_PER_WEEK * 86400
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
        return DurationObject(f"+{number}{unit}")

    @classmethod
    def from_expires_in_to_str(cls, expires_in: datetime) -> str:
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

    @classmethod
    def from_seconds(cls, seconds: int):
        duration = DurationObject(f"{seconds}s")
        return duration
