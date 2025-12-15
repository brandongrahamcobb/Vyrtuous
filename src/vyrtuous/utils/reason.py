''' reason.py The purpose of this program is to provide the Reason utility class.
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

class Reason:

    def __init__(self):
        self.action: Optional[str] = None
        self.action_types = ['append', 'overwrite', 'delete']
        self.old_reason: Optional[str] = None
        self.new_reason: Optional[str] = None
        
    def interpret_action(self, action_prefix: Optional[str] = None) -> Optional[str]:
        action = None
        if action_prefix == '+':
            action = self.action_types[0]
        if action_prefix == '=':
            action = self.action_types[1]
        if action_prefix == '-':
            action = self.action_types[2]
        return action    
        
    def load_action(self, action: Optional[str] = None) -> None:
        self.action = action
        
    def load_old_reason(self, old_reason: Optional[str] = 'No reason provided.') -> None:
        self.old_reason = old_reason
    
    def load_new_reason(self, new_reason: Optional[str] = 'No reason provided.') -> None:
        self.new_reason = new_reason
        
    def output_display(self) -> None:
        action = self.interpret_action(self.action)
        if action in self.action_types:
            match action:
                case 'append':
                    if self.new_reason:
                        new_reason = self.new_reason.strip()
                    else:
                        new_reason = ''
                    if not new_reason:
                        raise ValueError('No new reason to append was provided.')
                    if self.old_reason:
                        updated_reason = f'{self.old_reason}\n{self.new_reason}' 
                    else:
                        updated_reason = self.new_reason
                case 'overwrite':
                    updated_reason = self.new_reason.strip() if self.new_reason else ''
                    if not updated_reason:
                        raise ValueError('No new reason to overwrite was provided.')
                case 'delete':
                    updated_reason = None
                case _:
                    raise
        else:
            raise ValueError(f'Action does not match any of the valid types: {' '.join(self.action_types)}.')
        return updated_reason
