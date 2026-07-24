from dataclasses import dataclass, field

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
            "server-mutes": [],
            "summary": [],
            "survey": [],
            "streams": [],
            "text-mutes": [],
            "video-channels": [],
            "voice-mutes": [],
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


@dataclass
class PermissionGroup:
    name: str
    alias: str
    default: bool
    permissions: set[str] = field(default_factory=set)
    inheritance: list[str] = field(default_factory=list)


@dataclass
class PermissionNode:
    name: str
    path: str
    children: dict[str, "PermissionNode"] = field(default_factory=dict)
