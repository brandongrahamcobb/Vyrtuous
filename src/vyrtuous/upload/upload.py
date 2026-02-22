from dataclasses import dataclass


@dataclass(frozen=True)
class Upload:
    __tablename__ = "uploads"
    command_name: str
    arguments: str
    file_bytes: bytes
    filename: str
