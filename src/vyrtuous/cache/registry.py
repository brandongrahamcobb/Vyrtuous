import asyncio
import time
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from typing import Type, TypeVar, Union, cast

import discord
from cachetools import TTLCache

T = TypeVar("T")


class Registry:
    def __init__(self):
        self._services: dict[type, object] = {}

    def register(self, service: Union[object, tuple[object, ...]]) -> "Registry":
        if isinstance(service, tuple):
            for s in service:
                self._services[type(s)] = s
        else:
            self._services[type(service)] = service
        return self

    def get(self, t: Type[T]) -> T:
        return cast(T, self._services.get(t))


@dataclass
class MemberState:
    active: dict[int, tuple[str, datetime]] = field(default_factory=dict)
    flagged: defaultdict[int, set[int]] = field(
        default_factory=lambda: defaultdict(set)
    )
    invincible: defaultdict[int, set[int]] = field(
        default_factory=lambda: defaultdict(set)
    )


@dataclass
class ChannelState:
    deleted: set[str] = field(default_factory=set)
    join_log: defaultdict[tuple[int, int], list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    video: set[int] = field(default_factory=set)
    join_log_window: float = 300.0

    def should_warn(self, channel_id: int, member_id: int) -> bool:
        key = (channel_id, member_id)
        now = time.time()
        self.join_log[key] = [
            t for t in self.join_log[key] if now - t < self.join_log_window
        ]
        return len(self.join_log[key]) < 1

    def record(self, channel_id: int, member_id: int) -> None:
        self.join_log[(channel_id, member_id)].append(time.time())


@dataclass
class MessageHistoryState:
    cache: TTLCache = field(
        default_factory=lambda: TTLCache(maxsize=500, ttl=8 * 60 * 60)
    )
    reporters: defaultdict[int, set[int]] = field(
        default_factory=lambda: defaultdict(set)
    )


@dataclass
class VideoRoomState:
    tasks: dict[tuple[int, int], asyncio.Task] = field(default_factory=dict)
    cooldowns: dict[int, datetime] = field(default_factory=dict)

    def schedule(self, member: discord.Member, delay: int, coro) -> None:
        key = (member.guild.id, member.id)
        self.cancel(key)
        self.tasks[key] = asyncio.create_task(coro)

    def cancel(self, key: tuple[int, int]) -> None:
        task = self.tasks.pop(key, None)
        if task:
            task.cancel()

    def is_on_cooldown(self, member_id: int, cooldown: timedelta) -> bool:
        last = self.cooldowns.get(member_id)
        return last is not None and datetime.now(timezone.utc) - last < cooldown

    def set_cooldown(self, member_id: int) -> None:
        self.cooldowns[member_id] = datetime.now(timezone.utc)


@dataclass
class SystemResourcesState:
    cpu_seconds: list[float] = field(default_factory=list)
    rx_bytes: list[int] = field(default_factory=list)
    tx_bytes: list[int] = field(default_factory=list)
