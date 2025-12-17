''' time_to_complete.py The purpose of this program is to provide the TimeToComplete utility module.
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
from typing import Optional

class TimeToComplete:
 
    def __init__(self):
        pass

    def is_around_one_second(self, elapsed: float = 1.0):
        return (0.0 <= elapsed) and (elapsed <= 2.0)
            
    def time_elapsed_measurement(self, start: Optional[int], end: Optional[int]) -> float:
        if start > end:
            return abs(end - start)
        return end - start

