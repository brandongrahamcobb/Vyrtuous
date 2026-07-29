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

PERMISSION_TREE = {
    "other_guilds": [],
    "command": {
        "alias": ["create", "delete"],
        "utility": ["delete", "ping", "purge", "move"],
        "info": {
            "aliases": [],
            "automutes": [],
            "bans": [],
            "blacklists": [],
            "caps": [],
            "cogs": [],
            "debug": [],
            "flags": [],
            "guilds": [],
            "intents": [],
            "heroes": [],
            "overwrites": [],
            "roleid": [],
            "roles": [],
            "scope": ["channel", "guild", "member", "role"],
            "stats": [],
            "summary": [],
            "survey": [],
            "streams": [],
            "text-mutes": [],
            "vegans": [],
            "video-channels": [],
            "voice-mutes": ["auto", "click", "command", "server"],
        },
        "clear": {
            "category": [
                "alias",
                "automute",
                "ban",
                "flag",
                "stream",
                "text-mute",
                "video-channel",
                "voice_mute",
            ],
            "target": ["auto", "click", "command", "server"],
            "scope": ["all", "channel", "guild", "member", "role"],
        },
        "moderation": {
            "ban": [],
            "blacklist": [],
            "duration": [],
            "flag": [],
            "hero": [],
            "reason": [],
            "text-mute": [],
            "voice-mute": [
                "auto",
                "channel_mute",
                "command",
                "server",
            ],
            "unban": [],
            "unflag": [],
            "untext-mute": [],
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
            "load",
            "reload",
            "sync",
            "unload",
            "upload",
        ],
        "groups": ["grant", "revoke"],
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
