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
from typing import Optional

class Duration:
    
    now: datetime = datetime.now(timezone.utc)
            
    def __init__(self):
        self.base: Optional[datetime] = None
        self.number: Optional[int] = None
        self.prefix: Optional[str] = None
        self.suffix: Optional[str] = None
        self.day_suffix_types = ['d', 'day', 'days']
        self.hour_suffix_types = ['h', 'hr', 'hrs', 'hour', 'hours']
        self.minute_suffix_types = ['m', 'min', 'mins', 'minute', 'minutes']
        self.prefix_types = ['-', '+', '=']
        self.suffix_list = [self.day_suffix_types, self.hour_suffix_types, self.minute_suffix_types]

    def get_prefix(self, duration: Optional[str]) -> Optional[str]:
        if duration and duration[0] in self.prefix_types:
            prefix = duration[0]
        else:
            prefix = '+'
        return prefix
    
    def get_number(self, duration: Optional[str]) -> Optional[int]:
        number_str = ''
        if not duration:
            return 0 
        if duration and duration[0] in self.prefix_types:
            duration = duration[1:]
        for char in duration:
            if char.isdigit():
                number_str += char
            else:
                break
        if not number_str:
            raise ValueError('No numeric value found in duration string')
        number = int(number_str)
        return number
    
    def get_suffix(self, duration: Optional[str]) -> Optional[str]:
        if duration and duration[0] in self.prefix_types:
            duration = duration[1:]
        number_found = False
        suffix = ''
        for char in duration:
            if char.isdigit():
                number_found = True
            elif number_found:
                suffix += char
        if not suffix:
            suffix = 'h'
        suffix = suffix.lower()
        return suffix


    def interpret_days_suffix(self):
        if self.suffix in self.suffix_list[0]:
            return True
        else:
            return False
            
    def interpret_hours_suffix(self):
        if self.suffix in self.suffix_list[1]:
            return True
        else:
            return False
            
    def interpret_minutes_suffix(self):
        if self.suffix in self.suffix_list[2]:
            return True
        else:
            return False
    
    def interpret_sign(self) -> Optional[int]:
        if self.prefix == '-':
            return -1
        if self.prefix in ('+', '='): 
            return 1
        return 1
    
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
        self.load_prefix(prefix)
        number = self.get_number(duration_str)
        self.load_number(number)
        if self.validate_permanent():
            return
        suffix = self.get_suffix(duration_str)
        self.load_suffix(suffix)
            
    def load_number(self, number: Optional[int]) -> None:
        self.number = number
        
    def load_prefix(self, prefix: Optional[str]) -> None:
        self.prefix = prefix
        
    def load_suffix(self, suffix: Optional[str]) -> None:
        self.suffix = suffix
        
    def validate_base(self) -> None:
        if self.base is None:
           self.base = datetime.now(timezone.utc)
        
    def validate_number(self) -> None: 
        if self.number is None:
            raise ValueError('Duration.number is unset')
    
    def validate_permanent(self) -> bool:
        return self.number == 0
        
    def validate_prefix(self) -> None:
        if self.prefix is None:
            raise ValueError('Duration.prefix is unset')
        if self.prefix not in self.prefix_types:
            raise ValueError(f'Duration.prefix {self.prefix} is invalid')
    
    def validate_suffix(self) -> None:
        if self.suffix is None:
            raise ValueError('Duration.suffix is unset')
        if self.suffix not in sum(self.suffix_list, []):
            raise ValueError(f'Duration.suffix {self.suffix} is invalid')

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
        delta = self.build_timedelta()
        if self.base:
            base_time = self.base
        else:
            base_time = self.now
        return base_time + delta
        
    @classmethod
    def convert_timedelta_seconds(cls, seconds: int) -> str:
        if seconds == 0:
            return "0"  # or "permanent" depending on your convention
        sign = "+" if seconds > 0 else "-"
        seconds_abs = abs(seconds)
        days, remainder = divmod(seconds_abs, 86400)
        hours, remainder = divmod(remainder, 3600)
        minutes, _ = divmod(remainder, 60)
        if days > 0:
            return f"{sign}{days}d"
        elif hours > 0:
            return f"{sign}{hours}h"
        else:
            return f"{sign}{minutes}m"
        
    def output_datetime(self) -> Optional[datetime]:
        if self.validate_permanent():
            return None
        self.validate_base()
        self.validate_number()
        self.validate_prefix()
        self.validate_suffix()
        return self.compute_target_datetime()
    
    def output_display(self) -> Optional[str]:
        if self.validate_permanent():
            return "permanent"
        self.validate_prefix()
        self.validate_number()
        self.validate_suffix()
        sign = self.interpret_sign()
        number = self.number
        unit = self.interpret_unit()
        if sign > 0:
            return f"{number}{unit} from now"
        elif sign < 0:
            return f"{number}{unit} ago"
        else:
            return f"{number}{unit}"
    
    @classmethod
    def output_display_from_datetime(cls, dt: Optional[datetime]) -> Optional[str]:
        if dt is None:
            return "permanent"
        delta = dt - cls.now
        seconds = int(delta.total_seconds())
        if seconds == 0:
            return "now"
        sign = 1 if seconds > 0 else -1
        seconds_abs = abs(seconds)
        days = seconds_abs // 86400
        hours = seconds_abs // 3600
        minutes = seconds_abs // 60
        match True:
            case _ if days >= 365:
                number, unit = days // 365, "years"
            case _ if days >= 30:
                number, unit = days // 30, "months"
            case _ if days >= 7:
                number, unit = days // 7, "weeks"
            case _ if days >= 1:
                number, unit = days, "days"
            case _ if hours >= 1:
                number, unit = hours, "hours"
            case _:
                number, unit = minutes, "minutes"
        unit_text = unit.rstrip('s') if number == 1 else unit
        return f"{number}{unit_text} from now" if sign > 0 else f"{number}{unit_text} ago"
