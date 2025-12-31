''' duration.py The purpose of this program is to provide the Duration utility class.

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
from datetime import datetime, timedelta, timezone
from discord import app_commands
from discord.ext import commands
from typing import Optional
import discord

class DurationObject:
            
    DAY_UNITS = {'d', 'day', 'days'}
    HOUR_UNITS = {'h', 'hr', 'hrs', 'hour', 'hours'}
    MINUTE_UNITS = {'m', 'min', 'mins', 'minute', 'minutes'}
    SECOND_UNITS = {'s', 'sec', 'secs', 'second', 'seconds'}
    PREFIXES = {'+', '-', '='}

    UNIT_MAP = {**dict.fromkeys(DAY_UNITS, 'd'),
                **dict.fromkeys(HOUR_UNITS, 'h'),
                **dict.fromkeys(MINUTE_UNITS, 'm'),
                **dict.fromkeys(SECOND_UNITS, 's')}

    def __init__(self, duration: str):
        self._duration: str = ""
        self._prefix: Optional[str] = None
        self._number: Optional[int] = None
        self._unit: Optional[str] = None
        self._base: Optional[datetime] = None
        self.is_modification: bool = None
        self.duration = duration

    def __str__(self):
        if self.number == 0:
            return "permanent"
        return f"<t:{int(self.expires_at.timestamp())}:R>"

    @property
    def duration(self) -> str:
        return self._duration

    @duration.setter
    def duration(self, value: str):
        if value is None:
            raise commands.BadArgument("Duration string must be non-empty string")
        if isinstance(value, int):
            self._duration = str(value)
        else:
            self._duration = value.strip()
        self._parse()

    @property
    def prefix(self) -> str:
        return self._prefix or '+'

    @prefix.setter
    def prefix(self, value: str):
        if (value not in self.PREFIXES) or value is None:
            raise commands.BadArgument(f"Invalid prefix: {value}")
        self._prefix = value

    @property
    def number(self) -> int:
        if self._number is None:
            raise ValueError("Duration.number is unset")
        return self._number

    @number.setter
    def number(self, value: int):
        if not isinstance(value, int) or value < 0:
            raise commands.BadArgument("Duration.number must be a non-negative integer")
        self._number = value

    @property
    def unit(self) -> str:
        return self._unit or 'h'

    @unit.setter
    def unit(self, value: str):
        if value not in self.UNIT_MAP.values():
            raise commands.BadArgument(f"Invalid unit: {value}")
        self._unit = value

    @property
    def base(self) -> datetime:
        if self._base is None:
            self._base = datetime.now(timezone.utc)
        return self._base

    @base.setter
    def base(self, value: Optional[datetime]):
        self._base = value or datetime.now(timezone.utc)

    @property
    def sign(self) -> int:
        return -1 if self.prefix == '-' else 1

    @property
    def expires_at(self) -> datetime:
        return self.target_datetime()
    
    @classmethod
    def from_timedelta(cls, td: timedelta, prefix: str = '+') -> "DurationObject":
        total_seconds = int(td.total_seconds())
        if total_seconds % 86400 == 0:
            number, suffix = total_seconds // 86400, 'd'
        elif total_seconds % 3600 == 0:
            number, suffix = total_seconds // 3600, 'h'
        elif total_seconds % 60 == 0:
            number, suffix = total_seconds // 60, 'm'
        else:
            number, suffix = total_seconds, 's'
        obj = cls(f"{prefix}{number}{suffix}")
        return obj
    
    @classmethod
    def from_expires_at(cls, expires_at: datetime) -> "DurationObject":
        from datetime import datetime, timezone, timedelta
        now = datetime.now(timezone.utc)
        remaining = expires_at - now
        total_seconds = int(remaining.total_seconds())
        if total_seconds < 0:
            total_seconds = 0
        if total_seconds % 86400 == 0:
            number, unit = total_seconds // 86400, 'd'
        elif total_seconds % 3600 == 0:
            number, unit = total_seconds // 3600, 'h'
        elif total_seconds % 60 == 0:
            number, unit = total_seconds // 60, 'm'
        else:
            number, unit = total_seconds, 's'
        return cls(f"+{number}{unit}")
    
    @classmethod
    def from_seconds(cls, seconds: int):
        duration = DurationObject(f"{seconds}s")
        return duration
    
    def to_timedelta(self) -> timedelta:
        match self.unit:
            case 'd':
                return timedelta(days=self.number * self.sign)
            case 'h':
                return timedelta(hours=self.number * self.sign)
            case 'm':
                return timedelta(minutes=self.number * self.sign)
            case 's':
                return timedelta(seconds=self.number * self.sign)
            case _:
                raise ValueError(f"Unsupported unit: {self.unit}")

    def target_datetime(self, base: Optional[datetime] = None) -> datetime:
        base = base or self.base
        return base + self.to_timedelta()

    def discord_unix_timestamp(self, base: Optional[datetime] = None) -> int:
        return int(self.target_datetime(base).timestamp())

    def to_seconds(self) -> int:
        return int(self.to_timedelta().total_seconds())

    def _parse(self):
        s = self._duration.lower()
        if s == "0":
            self.number = 0
            self.unit = 'h'
            self.prefix = '+'
            return
        if s[0] in self.PREFIXES:
            self.prefix = s[0]
            self.is_modification = True
            s = s[1:]
        else:
            self.prefix = '+'
            self.is_modification = False
        num_str = ''
        for char in s:
            if char.isdigit():
                num_str += char
            else:
                break
        if not num_str:
            raise commands.BadArgument(f"No numeric duration found in '{self._duration}'")
        self.number = int(num_str)
        s = s[len(num_str):].strip()
        if not s:
            self.unit = 'h'
        else:
            unit = self.UNIT_MAP.get(s)
            if not unit:
                for known in self.UNIT_MAP.keys():
                    if s.startswith(known):
                        unit = self.UNIT_MAP[known]
                        break
            if not unit:
                raise commands.BadArgument(f"Invalid duration unit in '{self._duration}'")
            self.unit = unit


class Converter(commands.Converter):
    
    def __init__(self, duration_cls=DurationObject):
        self.duration_cls = duration_cls
    
    async def convert(self, ctx: commands.Context, arg):
        return self.duration_cls(arg).duration

class Transformer(app_commands.Transformer):

    def __init__(self, duration_cls=DurationObject):
        self.duration_cls = duration_cls

    async def transform(self, interaction: discord.Interaction, arg):
        return self.duration_cls(arg).duration
        
class Duration(Converter):
    def __init__(self):
        super().__init__(DurationObject)

class AppDuration(Transformer):
    def __init__(self):
        super().__init__(DurationObject)