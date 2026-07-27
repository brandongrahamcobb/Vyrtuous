import yaml
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import (
    PERMISSION_TREE,
    PermissionGroup,
    PermissionNode,
    PermissionScope,
)
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.permission_entry import PermissionEntry
from vyrtuous.inc.helpers import PATH_GROUPS
from vyrtuous.utils.users.moderator_service import HasEqualOrLowerRole

MODEL = PermissionEntry


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
            is_sysadmin=value.get("sysadmin", False),
            is_guild_owner=value.get("guild_owner", False),
            scope=PermissionScope(value["scope"]),
            permissions=set(value.get("permissions", [])),
            inheritance=[parent.lower() for parent in value.get("inheritance", [])],
        )
    return groups


def resolve_group_permissions(
    group_alias: str,
    groups: dict[str, PermissionGroup],
    visited: set[str] | None = None,
) -> set[str]:
    if visited is None:
        visited = set()
    if group_alias in visited:
        raise ValueError(f"Circular inheritance detected: {group_alias}")
    visited.add(group_alias)
    group = groups[group_alias]
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
    for group in groups:
        if has_group_permission(
            permission_state=permission_state,
            group_alias=group.alias,
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
    sysadmin_groups = set()
    guild_owner_groups = set()
    for group in permission_state.groups.values():
        if group.is_sysadmin:
            sysadmin_groups.add(group)
            if len(sysadmin_groups) > 1:
                group_names = [group.name for group in sysadmin_groups]
                raise ValueError(
                    f"Multiple groups ({', '.join(group_names)}) have sysadmin privileges. Only one group can have sysadmin permissions."
                )
        if group.is_guild_owner:
            guild_owner_groups.add(group)
            if len(guild_owner_groups) > 1:
                group_names = [group.name for group in guild_owner_groups]
                raise ValueError(
                    f"Multiple groups ({', '.join(group_names)}) have guild owner privileges. Only one group can have guild owner permissions."
                )
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
    group_alias: str,
    requested: str,
) -> bool:
    group = permission_state.groups.get(group_alias)
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
) -> set[PermissionGroup]:
    groups = set()
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    if is_sysadmin(member_snowflake=member_snowflake):
        for group in permission_state.groups.values():
            if group.is_sysadmin:
                groups.add(group)
    if guild_snowflake:
        for guild in bot.guilds:
            if guild.id == guild_snowflake and guild.owner_id == member_snowflake:
                for group in permission_state.groups.values():
                    if group.is_guild_owner:
                        groups.add(group)
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    entries = await database_factory.select(
        member_snowflake=member_snowflake,
        singular=False,
    )
    for entry in entries:
        if entry.guild_snowflake is not None:
            if entry.guild_snowflake != guild_snowflake:
                continue
        if entry.channel_snowflake is not None:
            if entry.channel_snowflake != channel_snowflake:
                continue
        if group := permission_state.groups.get(entry.group_alias, None) is not None:
            groups.add(group)
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
    for group_alias in groups:
        groups[group_alias].ancestors = resolve_ancestors(
            group_alias,
            groups,
        )


def resolve_ancestors(
    group_alias: str,
    groups: dict[str, PermissionGroup],
    visited: set[str] | None = None,
) -> set[str]:
    if visited is None:
        visited = set()
    if group_alias in visited:
        raise ValueError(f"Circular inheritance detected: {group_alias}")
    visited.add(group_alias)
    ancestors = set()
    for parent in groups[group_alias].inheritance:
        ancestors.add(parent)
        ancestors.update(
            resolve_ancestors(
                parent,
                groups,
                visited.copy(),
            )
        )
    return ancestors


def is_sysadmin(member_snowflake: int) -> bool:
    bot: DiscordBot = DiscordBot.get_instance()
    if int(bot.config["discord_owner_id"]) == member_snowflake:
        return True
    return False


async def resolve_all_assigned_groups(
    permission_state: PermissionState,
    member_snowflake: int,
) -> list[tuple[PermissionGroup, int | None, int | None]]:
    assigned: list[tuple[PermissionGroup, int | None, int | None]] = []
    bot: DiscordBot = DiscordBot.get_instance()
    if is_sysadmin(member_snowflake=member_snowflake):
        for group in permission_state.groups.values():
            if group.is_sysadmin:
                assigned.append((group, None, None))
    for guild in bot.guilds:
        if guild.owner_id == member_snowflake:
            for group in permission_state.groups.values():
                if group.is_guild_owner:
                    assigned.append((group, guild.id, None))
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    entries = await database_factory.select(
        member_snowflake=member_snowflake, singular=False
    )
    for entry in entries:
        group = next(
            (
                g
                for g in permission_state.groups.values()
                if g.alias == entry.group_alias
            ),
            None,
        )
        if group is None:
            continue
        assigned.append((group, entry.guild_snowflake, entry.channel_snowflake))
    return assigned


async def resolve_effective_group(
    permission_state: PermissionState,
    member_snowflake: int,
    channel_snowflake: int | None = None,
    guild_snowflake: int | None = None,
) -> PermissionGroup | None:
    bot: DiscordBot = DiscordBot.get_instance()
    if is_sysadmin(member_snowflake=member_snowflake):
        for group in permission_state.groups.values():
            if group.is_sysadmin:
                return group
    if guild_snowflake:
        for guild in bot.guilds:
            if guild.id == guild_snowflake and guild.owner_id == member_snowflake:
                for group in permission_state.groups.values():
                    if group.is_guild_owner:
                        return group
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    entries = await database_factory.select(
        member_snowflake=member_snowflake,
        singular=False,
    )
    global_group: PermissionGroup | None = None
    guild_group: PermissionGroup | None = None
    channel_group: PermissionGroup | None = None
    for entry in entries:
        group = permission_state.groups.get(entry.group_alias.lower())
        if group is None:
            continue
        if entry.channel_snowflake is not None:
            if entry.channel_snowflake == channel_snowflake:
                channel_group = group
                continue
        if entry.guild_snowflake is not None:
            if entry.guild_snowflake == guild_snowflake:
                guild_group = group
                continue
        if entry.channel_snowflake is None and entry.guild_snowflake is None:
            global_group = group
    return channel_group or guild_group or global_group


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
