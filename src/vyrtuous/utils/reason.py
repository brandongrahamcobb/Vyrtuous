
import inspect
from typing import Optional

class Reason:

    action_types = ['append', 'overwrite', 'delete']

    def __init__(self):
        self.action: Optional[str] = None
        self.old_reason: Optional[str] = None
        self.new_reason: Optional[str] = None
        
    # Load Methods
    def load_action(self, action: Optional[str] = None) -> None:
        self.action = action
        
    def load_old_reason(self, old_reason: Optional[str] = None) -> None:
        self.old_reason = old_reason
    
    def load_new_reason(self, new_reason: Optional[str] = None) -> None:
        self.new_reason = new_reason
        
    # Create Method
    def output_display(self) -> None:
        try:
            action = self.interpret_action(self.action)
            if action in self.action_types:
                match action:
                    case 'append':
                        new_reason = self.new_reason.strip() if self.new_reason else ''
                        if not new_reason:
                            raise ValueError('No new reason to append was provided.')
                        updated_reason = f'{self.old_reason}\n{self.new_reason}' if self.old_reason else self.new_reason
                    case 'overwrite':
                        updated_reason = self.new_reason.strip() if self.new_reason else ''
                        if not updated_reason:
                            raise ValueError('No new reason to overwrite was provided.')
                    case 'delete':
                        updated_reason = None
                    case _:
                        raise
            else:
                raise ValueError(f'Action does not match any of the valid types: {' '.join(action_types)}.')
        except Exception as e:
            func = inspect.currentframe().f_code.co_name
            raise RuntimeError(f'Unhandled exception in {func} -> {e}')
        return updated_reason
        
    # Interpreter Method
    def interpret_action(self, action_prefix: Optional[str] = None) -> Optional[str]:
        action = None
        if action_prefix == '+':
            action = self.action_types[0]
        if action_prefix == '=':
            action = self.action_types[1]
        if action_prefix == '-':
            action = self.action_types[2]
        return action
        
