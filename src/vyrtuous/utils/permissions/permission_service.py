import yaml

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import PERMISSION_TREE, PermissionGroup, PermissionNode
from vyrtuous.cache.registry import PermissionState
from vyrtuous.inc.helpers import PATH_GROUPS


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
            permissions=set(value.get("permissions", [])),
            inheritance=[parent.lower() for parent in value.get("inheritance", [])],
        )
    return groups


def resolve_group_permissions(
    group_name: str,
    groups: dict[str, PermissionGroup],
    visited: set[str] | None = None,
) -> set[str]:
    group_name = group_name.lower()
    if visited is None:
        visited = set()
    if group_name in visited:
        raise ValueError(f"Circular inheritance detected: {group_name}")
    visited.add(group_name)
    group = groups[group_name]
    permissions = set(group.permissions)
    for parent in group.inheritance:
        permissions.update(
            resolve_group_permissions(
                parent,
                groups,
                visited.copy(),
            )
        )
    return permissions


def validate_groups(
    permission_state: PermissionState,
) -> None:
    for group in permission_state.groups.values():
        for permission in group.permissions:
            node = permission.removeprefix("-")
            if node == "*":
                continue
            if node.endswith(".*"):
                node = node.removesuffix(".*")

            if not permission_state.exists(node):
                raise ValueError(
                    f'Group "{group.name}" has unknown permission "{permission}"'
                )

        for parent in group.inheritance:
            if parent not in permission_state.groups:
                raise ValueError(
                    f'Group "{group.name}" inherits unknown group "{parent}"'
                )


def populate():
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state = bot.registry.get(PermissionState)
    root, nodes = build_permission_tree(permissions=PERMISSION_TREE)
    groups = load_groups(path=PATH_GROUPS)
    permission_state.root = root
    permission_state.nodes = nodes
    permission_state.groups = groups
    validate_groups(permission_state=permission_state)
