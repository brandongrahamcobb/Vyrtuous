
from datetime import datetime, timedelta, timezone
from typing import List, Optional
import inspect

class Duration:

    day_suffix_types = ['d', 'day', 'days']
    hour_suffix_types = ['h', 'hr', 'hrs', 'hour', 'hours']
    minute_suffix_types = ['m', 'min', 'mins', 'minute', 'minutes']
    prefix_types = ['-', '+', '=']
    suffix_group = [day_suffix_types, hour_suffix_types, minute_suffix_types]
    permanent_duration_group = []
    for suffix_types in suffix_group:
        for suffix in suffix_types:
            permanent_duration_group.append('0' + suffix)
            
    def __init__(self):
        self.base: Optional[datetime] = None
        self.now: datetime = datetime.now(timezone.utc)
        self.number: Optional[int] = None
        self.prefix: Optional[str] = None
        self.suffix: Optional[str] = None
    
    # Loading Methods
    def load_number(self, number: Optional[int]) -> None:
        self.number = number
        pass
        
    def load_prefix(self, prefix: Optional[str]) -> None:
        self.prefix = prefix
        pass
        
    def load_suffix(self, suffix: Optional[str]) -> None:
        self.suffix = suffix
        pass
        
    def load_base(self, base: Optional[datetime]) -> None:
        self.base = base
        pass
    
    # Builder Methods
    def build_timedelta(self):
        sign = self.determine_sign()
        unit = self.determine_unit()
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
        base_time = self.base if self.base else self.now
        return base_time + delta
        
    def output_datetime(self) -> Optional[datetime]:
        try:
            if self.validate_permanent():
                return None
            self.validate_base()
            self.validate_number()
            self.validate_prefix()
            self.validate_suffix()
            return self.compute_target_datetime()
        except Exception as e:
            func = inspect.currentframe().f_code.co_name
            raise RuntimeError(f'Unhandled exception in {func} -> {e}')
    
    def output_display(self) -> Optional[str]:
        try:
            if self.validate_permanent():
                return "permanent"
            self.validate_prefix()
            self.validate_number()
            self.validate_suffix()
            sign = self.determine_sign()
            number = self.number
            unit = self.determine_unit()
            if number == 1:
                unit_text = unit.rstrip('s')
            else:
                unit_text = unit
            if sign > 0:
                return f"{number}{unit_text} from now"
            elif sign < 0:
                return f"{number}{unit_text} ago"
            else:
                return f"{number}{unit_text}"
        except Exception as e:
            func = inspect.currentframe().f_code.co_name
            raise RuntimeError(f'Unhandled exception in {func} -> {e}')
    
    def output_display_from_datetime(self, dt: Optional[datetime]) -> Optional[str]:
        if dt is None:
            return "permanent"
        now = datetime.now(timezone.utc)
        delta = dt - now
        seconds = int(delta.total_seconds())
        if seconds == 0:
            return "now"
        sign = 1 if seconds > 0 else -1
        seconds_abs = abs(seconds)
        days = seconds_abs // 86400
        hours = seconds_abs // 3600
        minutes = seconds_abs // 60
        if days >= 365:
            number = days // 365
            unit = "years"
        elif days >= 30:
            number = days // 30
            unit = "months"
        elif days >= 7:
            number = days // 7
            unit = "weeks"
        elif days >= 1:
            number = days
            unit = "days"
        elif hours >= 1:
            number = hours
            unit = "hours"
        else:
            number = minutes
            unit = "minutes"
        if number == 1:
            unit_text = unit.rstrip('s')
        else:
            unit_text = unit
        if sign > 0:
            return f"{number}{unit_text} from now"
        return f"{number}{unit_text} ago"

            
    # Determination Methods
    def determine_days_suffix(self):
        if self.suffix in self.suffix_group[0]:
            return True
        else:
            return False
            
    def determine_hours_suffix(self):
        if self.suffix in self.suffix_group[1]:
            return True
        else:
            return False
            
    def determine_minutes_suffix(self):
        if self.suffix in self.suffix_group[2]:
            return True
        else:
            return False
    
    def determine_sign(self) -> Optional[int]:
        if self.prefix == '-':
            return -1
        if self.prefix in ('+', '='): 
            return 1
        return 1
    
    def determine_unit(self) -> Optional[str]:
        if self.determine_days_suffix():
            return 'd'
        if self.determine_hours_suffix():
            return 'h'
        if self.determine_minutes_suffix():
           return 'm'
        raise ValueError(f'Invalid suffix {self.suffix}')
    
    # Getter Methods
    def get_prefix(self, duration: Optional[str]) -> Optional[str]:
        if duration and duration[0] in self.prefix_types:
            prefix = duration[0]
        else:
            prefix = '+'
        return prefix
    
    def get_number(self, duration: Optional[str]) -> Optional[int]:
        number_str = ''
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
        suffix = suffix.lower()
        return suffix
        
    def load_duration(self, duration: Optional[str]) -> None:
        prefix = self.get_prefix(duration)
        self.load_prefix(prefix)
        number = self.get_number(duration)
        self.load_number(number)
        suffix = self.get_suffix(duration)
        self.load_suffix(suffix)
        
    # Validation Methods
    def validate_base(self) -> None:
        if self.base is None:
           self.base = self.now
    
    def validate_permanent(self, duration: Optional[str] = None) -> bool:
        return duration in self.permanent_duration_group if duration else self.number == 0
        
    def validate_number(self) -> None: 
        if self.number is None:
            raise ValueError('Duration.number is unset')
            
    def validate_prefix(self) -> None:
        if self.prefix is None:
            raise ValueError('Duration.prefix is unset')
        if self.prefix not in self.prefix_types:
            raise ValueError(f'Duration.prefix {self.prefix} is invalid')
    
    def validate_suffix(self) -> None:
        if self.suffix is None:
            raise ValueError('Duration.suffix is unset')
        if self.suffix not in sum(self.suffix_group, []):
            raise ValueError(f'Duration.suffix {self.suffix} is invalid')
