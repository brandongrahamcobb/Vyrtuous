import yaml
from discord.ext import commands

from vyrtuous.cache.permissions import PERMISSION_TREE, PermissionGroup, PermissionNode
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.permission_level import PermissionLevel
from vyrtuous.inc.helpers import PATH_GROUPS
from vyrtuous.utils.users.moderator_service import HasEqualOrLowerRole

MODEL = PermissionLevel


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


def exists(permission_state: PermissionState, permission: str) -> bool:
    if permission.endswith(".*"):
        permission = permission[:-2]
    return permission in permission_state.nodes


async def has_permission(
    permission_state: PermissionState,
    member_snowflake: int,
    requested: str,
    channel_snowflake: int | None = None,
    guild_snowflake: int | None = None,
) -> bool:
    groups = await get_member_groups(
        member_snowflake=member_snowflake,
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
    )
    for group_name in groups:
        if has_group_permission(
            permission_state=permission_state,
            group_name=group_name,
            requested=requested,
        ):
            return True
    return False


async def has_permissions(
    permission_state: PermissionState,
    member_snowflake: int,
    requested: list[str],
    channel_snowflake: int | None = None,
    guild_snowflake: int | None = None,
) -> bool:
    for permission in requested:
        if not await has_permission(
            permission_state=permission_state,
            member_snowflake=member_snowflake,
            requested=permission,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        ):
            raise commands.CheckFailure("You do not have sufficient access to do that.")
    return True


def validate_groups(permission_state: PermissionState) -> None:
    for group in permission_state.groups.values():
        for permission in group.permissions:
            node = permission.removeprefix("-")
            if node == "*":
                continue
            if node.endswith(".*"):
                node = node.removesuffix(".*")
            if not exists(permission_state, node):
                raise ValueError(
                    f'Group "{group.name}" has unknown permission "{permission}"'
                )
        for parent in group.inheritance:
            if parent not in permission_state.groups:
                raise ValueError(
                    f'Group "{group.name}" inherits unknown group "{parent}"'
                )


def has_group_permission(
    permission_state: PermissionState,
    group_name: str,
    requested: str,
) -> bool:
    group = permission_state.groups.get(group_name.lower())
    if group is None:
        return False
    for permission in group.permissions:
        deny = permission.startswith("-")
        node = permission.removeprefix("-")
        if matches(node, requested):
            return not deny
    return False


def matches(
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


async def get_member_groups(
    member_snowflake: int,
    channel_snowflake: int | None = None,
    guild_snowflake: int | None = None,
) -> set[str]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    roles = await database_factory.select(
        member_snowflake=member_snowflake,
        singular=False,
    )
    groups = set()
    for role in roles:
        if role.guild_snowflake is not None:
            if role.guild_snowflake != guild_snowflake:
                continue
        if role.channel_snowflake is not None:
            if role.channel_snowflake != channel_snowflake:
                continue
        groups.add(role.level_name.lower())
    return groups


def populate(permission_state: PermissionState):
    root, nodes = build_permission_tree(permissions=PERMISSION_TREE)
    groups = load_groups(path=PATH_GROUPS)
    validate_groups(permission_state=permission_state)
    load_group_ancestors(groups)
    permission_state.root = root
    permission_state.nodes = nodes
    permission_state.groups = groups


def load_group_ancestors(
    groups: dict[str, PermissionGroup],
) -> None:
    for group_name in groups:
        groups[group_name].ancestors = resolve_ancestors(
            group_name,
            groups,
        )


def resolve_ancestors(
    group_name: str,
    groups: dict[str, PermissionGroup],
    visited: set[str] | None = None,
) -> set[str]:
    if visited is None:
        visited = set()
    if group_name in visited:
        raise ValueError(f"Circular inheritance detected: {group_name}")
    visited.add(group_name)
    ancestors = {group_name}
    for parent in groups[group_name].inheritance:
        ancestors.update(
            resolve_ancestors(
                parent,
                groups,
                visited.copy(),
            )
        )
    return ancestors


async def resolve_effective_group(
    permission_state: PermissionState,
    member_snowflake: int,
    channel_snowflake: int | None = None,
    guild_snowflake: int | None = None,
) -> PermissionGroup | None:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    roles = await database_factory.select(
        member_snowflake=member_snowflake,
        singular=False,
    )
    global_role: PermissionGroup | None = None
    guild_role: PermissionGroup | None = None
    channel_role: PermissionGroup | None = None
    for role in roles:
        group = permission_state.groups.get(role.level_name.lower())
        if group is None:
            continue
        if role.channel_snowflake is not None:
            if role.channel_snowflake == channel_snowflake:
                channel_role = group
                continue
        if role.guild_snowflake is not None:
            if role.guild_snowflake == guild_snowflake:
                guild_role = group
                continue
        if role.channel_snowflake is None and role.guild_snowflake is None:
            global_role = group
    return channel_role or guild_role or global_role


async def has_equal_or_lower_role(
    permission_state: PermissionState,
    author_snowflake: int,
    member_snowflake: int,
    channel_snowflake: int | None = None,
    guild_snowflake: int | None = None,
) -> None:
    author_group = await resolve_effective_group(
        permission_state=permission_state,
        member_snowflake=author_snowflake,
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
    )
    member_group = await resolve_effective_group(
        permission_state=permission_state,
        member_snowflake=member_snowflake,
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
    )
    if author_group is None or member_group is None:
        raise HasEqualOrLowerRole
    if (
        member_group.alias not in author_group.ancestors
        or author_group.alias == member_group.alias
    ):
        raise HasEqualOrLowerRole
