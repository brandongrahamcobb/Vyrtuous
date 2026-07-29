import yaml

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import (
    PERMISSION_TREE,
    PermissionGroup,
    PermissionNode,
    PermissionScope,
)
from vyrtuous.cache.registry import PermissionState
from vyrtuous.inc.helpers import PATH_GROUPS
from vyrtuous.permissions import validate_permissions
from vyrtuous.permissions import permission_service


def build_permission_tree(
    permissions: dict,
    parent: PermissionNode | None = None,
    path: tuple[str, ...] = (),
    nodes: dict[str, PermissionNode] | None = None,
) -> tuple[PermissionNode, dict[str, PermissionNode]]:
    if nodes is None:
        nodes = {}
    node = parent or PermissionNode(name="root", path="")
    for name, value in permissions.items():
        current_path = (*path, name)
        current = PermissionNode(
            name=name,
            path=".".join(current_path),
        )
        node.children[name] = current
        nodes[current.path] = current
        if isinstance(value, dict):
            build_permission_tree(
                value,
                current,
                current_path,
                nodes,
            )
        elif isinstance(value, list):
            for child in value:
                child_path = (*current_path, child)
                child_node = PermissionNode(
                    name=child,
                    path=".".join(child_path),
                )
                current.children[child] = child_node
                nodes[child_node.path] = child_node
    return node, nodes


def populate() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    root, nodes = build_permission_tree(permissions=PERMISSION_TREE)
    groups = load_groups(path=PATH_GROUPS)
    permission_state.root = root
    permission_state.nodes = nodes
    permission_state.groups = groups
    validate_permissions.validate_groups()
    load_group_ancestors(groups)


def load_groups(path: str) -> dict[str, PermissionGroup]:
    with open(path, "r") as file:
        data = yaml.safe_load(file)
    groups = {}
    for name, value in data["groups"].items():
        alias = value.get("alias", name.lower())
        groups[alias] = PermissionGroup(
            alias=alias,
            name=name,
            default=value.get("default", False),
            is_sysadmin=value.get("sysadmin", False),
            is_guild_owner=value.get("guild_owner", False),
            scope=PermissionScope(value["scope"]),
            permissions=set(value.get("permissions", [])),
            inheritance=[parent.lower() for parent in value.get("inheritance", [])],
        )
    return groups


def load_group_ancestors(
    groups: dict[str, PermissionGroup],
) -> None:
    for group_alias in groups:
        groups[group_alias].ancestors = permission_service.resolve_ancestors(
            group_alias,
            groups,
        )
