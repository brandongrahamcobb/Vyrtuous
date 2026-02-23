from dataclasses import dataclass, field
from datetime import datetime, timezone


@dataclass(frozen=True)
class Upload:
    __tablename__ = "uploads"
    command_name: str
    file_bytes: bytes
    filename: str
    tag: str
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
