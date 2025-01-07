import datetime
import os
import subprocess

def setup_backup_directory(backup_dir: str) -> str:
    """
    Ensure the backup directory exists, and return the path.
    """
    os.makedirs(backup_dir, exist_ok=True)
    return backup_dir

def perform_backup(db_user: str, db_name: str, db_host: str, backup_dir: str) -> str:
    """
    Perform a PostgreSQL database backup and return the backup file path.
    
    :param db_user: The PostgreSQL username.
    :param db_name: The name of the database to back up.
    :param db_host: The host where the database is running.
    :param backup_dir: The directory where backups are stored.
    :return: The path to the created backup file.
    """
    # Generate a timestamped filename
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    backup_file = os.path.join(backup_dir, f"backup_{timestamp}.sql")

    # Build the pg_dump command
    dump_command = [
        "pg_dump",
        "-U", db_user,
        "-h", db_host,
        "-d", db_name,
        "-F", "p",  # Plain-text format
        "-f", backup_file,
    ]

    # Run the backup command
    result = subprocess.run(
        dump_command,
        capture_output=True,
        text=True,
    )

    # Check if the backup was successful
    if result.returncode != 0:
        raise RuntimeError(f"Backup failed: {result.stderr}")

    return backup_file
