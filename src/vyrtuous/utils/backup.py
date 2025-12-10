
from datetime import datetime
from typing import Optional

import os
import subprocess

class Backup:

    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    def __init__(self):
        self.database: Optional[str] = os.getenv('POSTGRES_DB')
        self.directory: Optional[str] = '/app/backups'
        self.file_name: Optional[str] = os.path.join(self.directory, f'backup_{self.timestamp}.sql')
        self.host: Optional[str] = os.getenv('POSTGRES_HOST')
        self.password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
        self.user: Optional[str] = os.getenv('POSTGRES_USER')
        
    # Load Methods
    def load_database(self, database: Optional[str]):
        if database and databse.strip():
            self.database = database
        else:
            raise ValueError(f"Database cannot be {isinstance(database)}.")
        pass
        
    def load_directory(self, directory: Optional[str]):
        if directory and directory.strip():
             self.directory = directory
        else:
            raise ValueError(f"Directory cannot be {isinstance(directory)}.")
        pass
        
    def load_host(self, host: Optional[str]):
        if host and host.strip():
            self.host = host
        else:
            raise ValueError(f"Host cannot be {isinstance(host)}.")
        pass
        
    def load_password(self, password: Optional[str]):
        if password and password.strip():
            self.password = password
        else:
            raise ValueError(f"Password cannot be {isinstance(password)}.")
        pass
        
    def load_user(self, user: Optional[str]):
        if user and user.strip():
            self.user = user
        else:
            raise ValueError(f"User cannot be {isinstance(user)}.")
        pass
        
    # Directory Creation
    def create_backup_directory(self) -> None:
        os.makedirs(self.directory, exist_ok=True)
        return
        
    # Execute Backup
    def execute_backup(self) -> None:
        dump_command = [
            'pg_dump',
            '-U', self.user,
            '-h', self.host,
            '-d', self.database,
            '-F', 'p',
            '-f', self.file_name,
        ]
        env = os.environ.copy()
        env['PGPASSWORD'] = self.password
        result = subprocess.run(
            dump_command,
            capture_output=True,
            text=True,
            env=env,
        )
        if result.returncode != 0:
            raise RuntimeError(f'Backup failed: {result.stderr}')
        return result

