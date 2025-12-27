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
            
    day_suffix_types = ['d', 'day', 'days']
    hour_suffix_types = ['h', 'hr', 'hrs', 'hour', 'hours']
    minute_suffix_types = ['m', 'min', 'mins', 'minute', 'minutes']
    prefix_types = ['-', '+', '=']
    suffix_list = [day_suffix_types, hour_suffix_types, minute_suffix_types]

    def __init__(self, duration_str: str):
        self.duration_str = duration_str
        self.base: Optional[datetime] = None
        self.number: Optional[int] = None
        self.prefix: Optional[str] = None
        self.suffix: Optional[str] = None

    def __str__(self):
        return self._duration_str

    @property
    def duration_str(self) -> str:
        return self._duration_str

    @duration_str.setter
    def duration_str(self, duration: str):
        if not isinstance(duration, str) or not duration:
            raise commands.BadArgument("Duration string is invalid.")
        self._duration_str = duration.strip()
        self._parse_duration()

    @property
    def prefix(self) -> str:
        if self._prefix is None:
            return '+'
        if self._prefix not in self.prefix_types:
            raise ValueError(f"Duration.prefix {self._prefix} is invalid")
        return self._prefix

    @prefix.setter
    def prefix(self, value: str):
        if value not in self.prefix_types:
            raise commands.BadArgument(f"Duration.prefix {value} is invalid")
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
    def suffix(self) -> str:
        if self._suffix is None:
            return 'h'
        if self._suffix not in sum(self.suffix_list, []):
            raise ValueError(f"Duration.suffix {self._suffix} is invalid")
        return self._suffix

    @suffix.setter
    def suffix(self, value: str):
        if value.lower() not in sum(self.suffix_list, []):
            raise commands.BadArgument(f"Duration.suffix {value} is invalid")
        self._suffix = value.lower()

    @property
    def base(self) -> datetime:
        if self._base is None:
            self._base = datetime.now(timezone.utc)
        return self._base

    @base.setter
    def base(self, value: Optional[datetime]):
        self._base = value or datetime.now(timezone.utc)
    
    def interpret_sign(self) -> Optional[int]:
        return -1 if self.prefix == '-' else 1
    
    def interpret_unit(self) -> Optional[str]:
        if self.interpret_days_suffix():
            return 'd'
        if self.interpret_hours_suffix():
            return 'h'
        if self.interpret_minutes_suffix():
           return 'm'
        raise ValueError(f'Invalid suffix {self.suffix}')
        
    def load_base(self, base: Optional[datetime]) -> None:
        self.base = base
        pass
    
    def load_from_combined_duration_str(self, duration_str: Optional[str]) -> None:
        prefix = self.get_prefix(duration_str)
        self.prefix = prefix
        number = self.get_number(duration_str)
        self.number = number
        if self.validate_permanent():
            return
        suffix = self.get_suffix(duration_str)
        self.suffix = suffix
            
    def build_timedelta(self):
        sign = self.interpret_sign()
        unit = self.interpret_unit()
        match unit:
            case 'd':
                return timedelta(days=self.number * sign)
            case 'h':
                return timedelta(hours=self.number * sign)
            case 'm':
                return timedelta(minutes=self.number * sign)
            case _:
                return None
                
    def compute_target_datetime(self) -> Optional[datetime]:
        base = base or self.base
        return base + self.build_timedelta()
    
    def output_display(self) -> str:
        if self.number == 0:
            return "permanent"
        sign = self.interpret_sign()
        unit = self.interpret_unit()
        number = self.number
        if sign > 0:
            return f"{number}{unit} from now"
        else:
            return f"{number}{unit} ago"

class Converter(commands.Converter):
    
    def __init__(self, duration=DurationObject):
        self.duration = duration
    
    async def convert(self, ctx: commands.Context, arg):
        return self.duration(arg).duration

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