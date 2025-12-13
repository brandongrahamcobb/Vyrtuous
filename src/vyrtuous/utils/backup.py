
from datetime import datetime
from typing import Optional

import os
import subprocess

class Backup:

    def __init__(self, directory: Optional[str]):
        self.database: Optional[str] = os.getenv('POSTGRES_DB')
        self.directory = directory
        self.host: Optional[str] = os.getenv('POSTGRES_HOST')
        self.password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
        self.timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.file_name: Optional[str] = os.path.join(self.directory, f'backup_{self.timestamp}.sql')
        self.user: Optional[str] = os.getenv('POSTGRES_USER')
        
    def create_backup_directory(self) -> None:
        os.makedirs(self.directory, exist_ok=True)
        return
        
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

