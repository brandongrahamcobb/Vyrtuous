from dataclasses import dataclass, field

from discord.ext import commands

PERMISSION_TREE = {
    "list": {
        "aliases": {
            "scope": ["channel", "guild"],
            "other_guilds": [],
        },
        "automutes": {
            "scope": ["channel", "guild"],
            "other_guilds": [],
        },
        "bans": {
            "scope": ["channel", "guild", "member"],
            "other_guilds": [],
        },
        "blacklists": {
            "scope": ["channel", "guild", "member"],
            "other_guilds": [],
        },
        "caps": {"scope": ["channel", "guild"], "other_guilds": []},
        "flags": {
            "scope": ["channel", "guild", "member"],
            "other_guilds": [],
        },
        "heroes": {"scope": ["guild", "member"], "other_guilds": []},
        "groups": {"scope": ["group", "guild", "member"], "other_guilds": []},
        "roles": {"scope": ["role"], "other_guilds": []},
        "server-mutes": {
            "scope": ["guild", "member"],
            "other_guilds": [],
        },
        "summary": {"scope": ["member"], "other_guilds": []},
        "survey": {"scope": ["channel"], "other_guilds": []},
        "streams": {"scope": ["channel", "guild"], "other_guilds": []},
        "text-mutes": {
            "scope": ["channel", "guild", "member"],
            "other_guilds": [],
        },
        "video-channels": {
            "scope": ["channel", "guild"],
            "other_guilds": [],
        },
        "voice-mutes": {
            "scope": ["channel", "guild", "member"],
            "other_guilds": [],
        },
    },
    "clear": {
        "category": {
            "automute": ["channel", "guild"],
            "scope": {
                "channel": [],
                "guild": [],
            },
        },
        "ban": {
            "scope": {
                "channel": [],
                "guild": [],
                "member": [],
            }
        },
        "flag": {
            "scope": {
                "channel": [],
                "guild": [],
                "member": [],
            }
        },
        "text-mute": {
            "scope": {
                "channel": [],
                "guild": [],
                "member": [],
            }
        },
        "voice_mute": {
            "target": ["auto", "click", "command", "server"],
            "scope": {
                "channel": [],
                "guild": [],
                "member": [],
            },
        },
        "group": [],
    },
    "movement": ["channel_move"],
    "info": {
        "overwrites": {"scope": ["channel"], "other_guilds": []},
        "permissions": {"scope": ["channel"], "other_guilds": []},
        "roleid": {"scope": ["role"], "other_guilds": []},
    },
    "moderation": {
        "ban": ["capped", "uncapped"],
        "blacklist": [],
        "unban": [],
        "duration": ["uncapped"],
        "flag": [],
        "unflag": [],
        "hero": [],
        "reason": [],
        "text-mute": ["uncapped"],
        "untext-mute": [],
        "voice-mute": {
            "auto": ["uncapped"],
            "channel_mute": ["uncapped"],
            "command": ["uncapped"],
            "server": [],
        },
        "unvoice-mute": ["channel_unmute"],
    },
    "channel": {
        "automute": [],
        "cap": [],
        "stream": [],
        "video": [],
    },
    "dev": {
        "backup": [],
        "cogs": [],
        "debug": [],
        "load": [],
        "ping": [],
        "reload": [],
        "sync": [],
        "stats": [],
        "unload": [],
        "upload": [],
    },
}


@dataclass
class PermissionGroup:
    name: str
    permissions: set[str]


@dataclass
class PermissionNode:
    name: str
    path: str
    children: dict[str, "PermissionNode"] = field(default_factory=dict)


def build_permission_tree(
    permissions: dict,
    parent: PermissionNode | None = None,
    path: tuple[str, ...] = (),
) -> PermissionNode:
    node = parent or PermissionNode(name="root", path="")
    for name, value in permissions.items():
        current_path = (*path, name)
        current = PermissionNode(
            name=name,
            path=".".join(current_path),
        )
        node.children[name] = current
        if isinstance(value, dict):
            build_permission_tree(
                value,
                current,
                current_path,
            )
        elif isinstance(value, list):
            for child in value:
                child_path = (*current_path, child)
                current.children[child] = PermissionNode(
                    name=child,
                    path=".".join(child_path),
                )
    return node


class PermissionObject:

    def __init__(self, permission: str):
        self.__permission = permission

    @property
    def permission(self) -> str:
        return self.__permission


class Converter(commands.Converter):

    def __init__(self, permission_cls=PermissionObject):
        self.permission_cls = permission_cls
        self.__identifiers = [obj.identifier for obj in CLASSES]

    async def convert(self, ctx: commands.Context, argument: str) -> PermissionObject:
        categories = self.__identifiers
        for extra in EXTRA_CATEGORIES:
            categories.append(extra)
        if argument not in categories:
            raise commands.CheckFailure(f"Invalid category specified: ({argument}).")
        return self.permission_cls(argument)


class Transformer(app_commands.Transformer):

    def __init__(self, permission_cls=PermissionObject):
        self.permission_cls = permission_cls
        self.__identifiers = [obj.identifier for obj in CLASSES]

    async def transform(
        self, interaction: discord.Interaction, arg: str
    ) -> PermissionObject:
        categories = self.__identifiers
        for extra in EXTRA_CATEGORIES:
            categories.append(extra)
        if arg not in categories:
            raise app_commands.CheckFailure(f"Invalid category specified: ({arg}).")
        return self.permission_cls(arg)


class Permission(Converter):
    def __init__(self):
        super().__init__(PermissionObject)


class AppPermission(Transformer):
    def __init__(self):
        super().__init__(PermissionObject)
