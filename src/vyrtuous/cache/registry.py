import asyncio
import time
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from typing import Self, Type, TypeVar, Union, cast

import discord
from cachetools import TTLCache

from vyrtuous.cache.permissions import PermissionGroup, PermissionNode

T = TypeVar("T")


class Registry:
    def __init__(self) -> None:
        self._services: dict[type, object] = {}

    def register(self, service: Union[object, tuple[object, ...]]) -> Self:
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
    automuted: dict[int, set[int]] = field(default_factory=lambda: defaultdict(set))
    flagged: defaultdict[int, set[int]] = field(
        default_factory=lambda: defaultdict(set)
    )
    invincible: defaultdict[int, set[int]] = field(
        default_factory=lambda: defaultdict(set)
    )


@dataclass
class ChannelState:
    joined_at: dict[tuple[int, int, int], float] = field(default_factory=dict)
    join_log: defaultdict[tuple[int, int, int], list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    video: set[int] = field(default_factory=set)

    def should_notify(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        timeout: float,
    ) -> bool:
        key = (channel_snowflake, guild_snowflake, member_snowflake)
        now = time.time()
        self.join_log[key] = [t for t in self.join_log[key] if now - t < timeout]
        return len(self.join_log[key]) < 1

    def record(
        self, channel_snowflake: int, guild_snowflake: int, member_snowflake: int
    ) -> None:
        self.join_log[(channel_snowflake, guild_snowflake, member_snowflake)].append(
            time.time()
        )


@dataclass
class MessageHistoryState:
    cache: TTLCache = field(
        default_factory=lambda: TTLCache(maxsize=500, ttl=8 * 60 * 60)
    )
    reporters: defaultdict[int, set[int]] = field(
        default_factory=lambda: defaultdict(set)
    )


@dataclass
class VideoChannelState:
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


@dataclass
class PermissionState:
    root: PermissionNode | None = None
    nodes: dict[str, PermissionNode] = field(default_factory=dict)
    groups: dict[str, PermissionGroup] = field(default_factory=dict)

    def exists(self, permission: str) -> bool:
        if permission.endswith(".*"):
            permission = permission[:-2]
        return permission in self.nodes

    def has_permission(
        self,
        group_name: str,
        requested: str,
    ) -> bool:
        group = self.groups.get(group_name.lower())
        if group is None:
            return False
        for permission in group.permissions:
            deny = permission.startswith("-")
            node = permission.removeprefix("-")
            if self.matches(node, requested):
                return not deny
        return False

    def matches(
        self,
        granted: str,
        requested: str,
    ) -> bool:
        granted_parts = granted.split(".")
        requested_parts = requested.split(".")
        if len(granted_parts) > len(requested_parts):
            return False
        for index, part in enumerate(granted_parts):
            if part == "*":
                return True
            if part != requested_parts[index]:
                return False
        return len(granted_parts) == len(requested_parts)
