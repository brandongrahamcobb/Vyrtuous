"""!/bin/python3
permissions.py The purpose of this program is to provide the permissions model.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from dataclasses import dataclass, field
from enum import Enum

# "groups": {"scope": ["group", "guild", "member"]},

PERMISSION_TREE = {
    "other_guilds": [],
    "uncapped": [],
    "command": {
        "listing": {
            "aliases": [],
            "automutes": [],
            "bans": [],
            "blacklists": [],
            "caps": [],
            "flags": [],
            "heroes": [],
            "roles": [],
            "scope": ["channel", "guild", "member", "role"],
            "summary": [],
            "survey": [],
            "streams": [],
            "text-mutes": [],
            "video-channels": [],
            "voice-mutes": ["auto", "click", "command", "server"],
        },
        "clear": {
            "category": ["automute", "ban", "flag", "text-mute", "voice_mute"],
            "target": ["auto", "click", "command", "server"],
            "scope": ["all", "channel", "guild", "member", "role"],
        },
        "movement": ["channel_move"],
        "info": [
            "overwrites",
            "intents",
            "roleid",
            "ping",
        ],
        "moderation": {
            "ban": [],
            "blacklist": [],
            "unban": [],
            "duration": [],
            "flag": [],
            "unflag": [],
            "hero": [],
            "reason": [],
            "text-mute": [],
            "untext-mute": [],
            "voice-mute": [
                "auto",
                "channel_mute",
                "command",
                "server",
            ],
            "unvoice-mute": ["channel_unmute"],
            "uncapped": [],
        },
        "channel": [
            "automute",
            "cap",
            "stream",
            "video",
        ],
        "dev": [
            "backup",
            "cogs",
            "debug",
            "load",
            "reload",
            "sync",
            "stats",
            "unload",
            "upload",
        ],
    },
}


class PermissionScope(str, Enum):
    GLOBAL = "global"
    GUILD = "guild"
    CHANNEL = "channel"


@dataclass
class PermissionGroup:
    name: str
    alias: str
    scope: PermissionScope
    is_guild_owner: bool = False
    is_sysadmin: bool = False
    default: bool = False
    permissions: set[str] = field(default_factory=set)
    inheritance: list[str] = field(default_factory=list)
    ancestors: set[str] = field(default_factory=set)


@dataclass
class PermissionNode:
    name: str
    path: str
    children: dict[str, "PermissionNode"] = field(default_factory=dict)
