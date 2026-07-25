"""!/bin/python3
duration.py The purpose of this program is to provide the Duration properties class.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from datetime import datetime, timedelta, timezone
from typing import Self

import discord
from discord import app_commands
from discord.ext import commands

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
UNIT_ORDER = [
    ("y", 86400 * DAYS_PER_YEAR),
    ("w", 86400 * DAYS_PER_WEEK),
    ("d", 86400),
    ("h", 3600),
    ("m", 60),
    ("s", 1),
]
UNIT_SECONDS = {
    "s": 1,
    "m": 60,
    "h": 3600,
    "d": 86400,
    "w": 604800,
    "y": 31536000,
}


class DurationObject:
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


def _largest_unit(total_seconds: int):
    for unit, seconds in UNIT_ORDER:
        if total_seconds % seconds == 0:
            return total_seconds // seconds, unit
    return total_seconds, "s"


class DurationBuilder:
    def __init__(self):
        self.__duration = DurationObject(number=8, prefix="", sign=1, unit="h")

    def load(self, duration: DurationObject) -> Self:
        self.__duration = duration
        return self

    def parse(self, value) -> Self:
        if not value:
            self.__duration = DurationObject(number=0, unit="", prefix="", sign=1)
            return self
        s = value.lower().strip()
        if s == "0":
            number = 0
            unit = "h"
            prefix = ""
            sign = 1
            self.__duration = DurationObject(
                number=number, prefix=prefix, sign=sign, unit=unit
            )
            return self
        if s[0] in "+-":
            sign = 1 if s[0] == "+" else -1
            s = s[1:]
        else:
            sign = 1
        if s and s[0] in PREFIXES:
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
            value = UNIT_MAP.get(s, None)
            if not value:
                for known in UNIT_MAP.keys():
                    if s.startswith(known):
                        unit = UNIT_MAP[known]
                        break
                raise commands.BadArgument(f"Invalid duration unit in '{value}'")
            else:
                unit = value
        self.__duration = DurationObject(
            number=number, unit=unit, prefix=prefix, sign=sign
        )
        return self

    def from_seconds(self, seconds) -> Self:
        self.parse(f"{seconds}s")
        return self

    def from_timedelta(self, td, prefix="+") -> Self:
        total_seconds = int(td.total_seconds())
        number, unit = _largest_unit(total_seconds)
        self.parse(f"{prefix}{number}{unit}")
        return self

    def to_seconds(self) -> int:
        unit_seconds = UNIT_SECONDS.get(getattr(self.__duration, "unit", "h"), 3600)
        return (
            getattr(self.__duration, "sign", 1)
            * getattr(self.__duration, "number", 0)
            * unit_seconds
        )

    def to_timedelta(self) -> timedelta:
        match self.__duration.unit:
            case "y":
                return timedelta(
                    days=self.__duration.number * DAYS_PER_YEAR * self.__duration.sign
                )
            case "w":
                return timedelta(
                    days=self.__duration.number * DAYS_PER_WEEK * self.__duration.sign
                )
            case "d":
                return timedelta(days=self.__duration.number * self.__duration.sign)
            case "h":
                return timedelta(hours=self.__duration.number * self.__duration.sign)
            case "m":
                return timedelta(minutes=self.__duration.number * self.__duration.sign)
            case "s":
                return timedelta(seconds=self.__duration.number * self.__duration.sign)
            case _:
                raise ValueError(f"Unsupported unit: {self.__duration.unit}")

    def to_expires_in(self, base: datetime | None = None) -> datetime | None:
        if not self.__duration.number:
            return None
        base = base or datetime.now(timezone.utc)
        return base + self.to_timedelta()

    def from_timestamp(self, expires_in: datetime) -> Self:
        now = datetime.now(timezone.utc)
        remaining = expires_in - now
        total_seconds = max(0, int(remaining.total_seconds()))
        number, unit = _largest_unit(total_seconds)
        self.parse(f"+{number}{unit}")
        return self

    def as_str(self) -> str:
        return f"{self.__duration.prefix}{self.__duration.number}{self.__duration.unit}"

    def build(self) -> DurationObject:
        return self.__duration

    def to_unix_ts(self, base: datetime | None = None) -> str:
        expires_in = self.to_expires_in(base)
        if expires_in is None or self.__duration.number == 0:
            return "permanent"
        return f"<t:{int(expires_in.timestamp())}:R>"


class DurationWrapper:

    def __init__(self, duration: DurationObject):
        self.__duration = duration

    @property
    def duration(self) -> DurationObject:
        return self.__duration


class Converter(commands.Converter):

    def __init__(self, duration_cls=DurationWrapper):
        self.duration_cls = duration_cls

    async def convert(self, ctx: commands.Context, argument: str) -> DurationWrapper:
        duration_builder = DurationBuilder()
        try:
            duration = duration_builder.parse(argument).build()
            return self.duration_cls(duration)
        except ValueError:
            raise commands.CheckFailure("This command must input a valid duration.")


class Transformer(app_commands.Transformer):

    def __init__(self, duration_cls=DurationWrapper):
        self.duration_cls = duration_cls

    async def transform(
        self, interaction: discord.Interaction, arg: str
    ) -> DurationWrapper:
        duration_builder = DurationBuilder()
        try:
            duration = duration_builder.parse(arg).build()
            return self.duration_cls(duration)
        except ValueError:
            raise app_commands.CheckFailure("This command must input a valid duration.")


class Duration(Converter):
    def __init__(self):
        super().__init__(DurationWrapper)


class AppDuration(Transformer):
    def __init__(self):
        super().__init__(DurationWrapper)
