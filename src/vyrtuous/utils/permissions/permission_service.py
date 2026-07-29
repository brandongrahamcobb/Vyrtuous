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
from vyrtuous.utils.errors.error import HasEqualOrLowerRole, NotGuildOwner

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


async def any_group_has_permission(
    permission_state: PermissionState,
    member_snowflake: int,
    requested: str,
) -> bool:
    groups_list = []
    bot: DiscordBot = DiscordBot.get_instance()
    for guild in bot.guilds:
        for channel in guild.channels:
            groups = await get_member_groups(
                member_snowflake=member_snowflake,
                channel_snowflake=channel.id,
                guild_snowflake=guild.id,
            )
            for group in groups:
                if group not in groups_list:
                    groups_list.append(group)
    for group in groups_list:
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


async def any_group_has_permissions(
    permission_state: PermissionState,
    member_snowflake: int,
    requested: list[str],
) -> bool:
    for permission in requested:
        if not await any_group_has_permission(
            permission_state=permission_state,
            member_snowflake=member_snowflake,
            requested=permission,
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
) -> list[PermissionGroup]:
    groups = []
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    if is_sysadmin(member_snowflake=member_snowflake):
        for group in permission_state.groups.values():
            if group.is_sysadmin:
                if group not in groups:
                    groups.append(group)
    if guild_snowflake:
        for guild in bot.guilds:
            if guild.id == guild_snowflake and guild.owner_id == member_snowflake:
                for group in permission_state.groups.values():
                    if group.is_guild_owner:
                        if group not in groups:
                            groups.append(group)
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
            if group not in groups:
                groups.append(group)
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
        return None
    if (
        member_group.alias not in author_group.ancestors
        or author_group.alias == member_group.alias
    ):
        raise HasEqualOrLowerRole(target_rank=member_group.alias)


async def is_guild_owner(
    member_snowflake: int, guild_snowflake: int | None = None
) -> bool:
    bot: DiscordBot = DiscordBot.get_instance()
    if guild_snowflake is None:
        for guild in bot.guilds:
            if guild and guild.owner_id == member_snowflake:
                return True
    else:
        guild = bot.get_guild(guild_snowflake)
        if guild and guild.owner_id == member_snowflake:
            return True
    raise NotGuildOwner(guild_snowflake=guild_snowflake)


# async def survey(channel_snowflake: int, guild_snowflake: int) -> list[discord.Embed]:
#     bot: DiscordBot = DiscordBot.get_instance()
#     guild = bot.get_guild(guild_snowflake)
#     if guild is None:
#         raise GuildNotFound(str(guild_snowflake))
#     channel = guild.get_channel(channel_snowflake)
#     if channel is None:
#         raise ChannelNotFound(str(channel_snowflake))
#     if not isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
#         raise commands.CheckFailure("This command must target a valid channel.")
#     chunk_size, pages = 7, []
#     (
#         sysadmins,
#         developers,
#         guild_owners,
#         administrators,
#         coordinators,
#         moderators,
#     ) = ([], [], [], [], [], [])
#
#     for member in channel.members:
#         try:
#             if await sysadmin_service.is_sysadmin(member_snowflake=member.id):
#                 sysadmins.append(member)
#         except commands.CheckFailure as e:
#             bot.logger.warning(str(e).capitalize())
#         try:
#             if await developer_service.is_developer(member_snowflake=member.id):
#                 developers.append(member)
#         except commands.CheckFailure as e:
#             bot.logger.warning(str(e).capitalize())
#         try:
#             if await guild_owner_service.is_guild_owner(
#                 guild_snowflake=channel.guild.id, member_snowflake=member.id
#             ):
#                 guild_owners.append(member)
#         except commands.CheckFailure as e:
#             bot.logger.warning(str(e).capitalize())
#         try:
#             if await administrator_service.is_administrator(
#                 guild_snowflake=channel.guild.id, member_snowflake=member.id
#             ):
#                 administrators.append(member)
#         except commands.CheckFailure as e:
#             bot.logger.warning(str(e).capitalize())
#         try:
#             if await coordinator_service.is_coordinator(
#                 channel_snowflake=channel.id,
#                 guild_snowflake=channel.guild.id,
#                 member_snowflake=member.id,
#             ):
#                 coordinators.append(member)
#         except commands.CheckFailure as e:
#             bot.logger.warning(str(e).capitalize())
#         try:
#             if await is_moderator(
#                 channel_snowflake=channel.id,
#                 guild_snowflake=channel.guild.id,
#                 member_snowflake=member.id,
#             ):
#                 moderators.append(member)
#         except commands.CheckFailure as e:
#             bot.logger.warning(str(e).capitalize())
#     sysadmins_chunks = [
#         sysadmins[i : i + chunk_size] for i in range(0, len(sysadmins), chunk_size)
#     ]
#     guild_owners_chunks = [
#         guild_owners[i : i + chunk_size]
#         for i in range(0, len(guild_owners), chunk_size)
#     ]
#     developers_chunks = [
#         developers[i : i + chunk_size] for i in range(0, len(developers), chunk_size)
#     ]
#     administrators_chunks = [
#         administrators[i : i + chunk_size]
#         for i in range(0, len(administrators), chunk_size)
#     ]
#     coordinators_chunks = [
#         coordinators[i : i + chunk_size]
#         for i in range(0, len(coordinators), chunk_size)
#     ]
#     moderators_chunks = [
#         moderators[i : i + chunk_size] for i in range(0, len(moderators), chunk_size)
#     ]
#     roles_chunks = [
#         ("Sysadmins", sysadmins, sysadmins_chunks),
#         ("Developers", developers, developers_chunks),
#         ("Guild Owners", guild_owners, guild_owners_chunks),
#         ("Administrators", administrators, administrators_chunks),
#         ("Coordinators", coordinators, coordinators_chunks),
#         ("Moderators", moderators, moderators_chunks),
#     ]
#     max_pages = max(len(c[2]) for c in roles_chunks)
#     for page in range(max_pages):
#         embed = discord.Embed(
#             title=f"{emojis.get_random_emoji()} Survey results for {channel.name}",
#             description=f"Total surveyed: {len(channel.members)}",
#             color=discord.Color.blurple(),
#         )
#         for role_name, role_list, chunks in roles_chunks:
#             chunk = chunks[page] if page < len(chunks) else []
#             embed.add_field(
#                 name=f"{role_name} ({len(chunk)}/{len(role_list)})",
#                 value=", ".join(u.mention for u in chunk) if chunk else "*None*",
#                 inline=False,
#             )
#         pages.append(embed)
#     return pages
#
# async def log_mod(
#     author_snowflake: int,
#     channel_snowflake: int,
#     display: bool,
#     guild_snowflake: int,
#     member_snowflake: int,
#     message_snowflake: int | None,
#     message_channel_snowflake: int | None,
# ):
#     duration = DurationObject(number=0, prefix="", sign=1, unit="")
#     is_channel_scope = None
#     reason = None
#     role_snowflake = None
#     target = None
#     await data_builder.save_data(
#         author_snowflake=author_snowflake,
#         channel_snowflake=channel_snowflake,
#         duration=duration,
#         guild_snowflake=guild_snowflake,
#         identifier="mod",
#         member_snowflake=member_snowflake,
#         reason=reason or "No reason provided.",
#         role_snowflake=role_snowflake or None,
#         target=target or None,
#     )
#     if display:
#         await stream_service.send_log(
#             author_snowflake=author_snowflake,
#             channel_snowflake=channel_snowflake,
#             identifier="mod",
#             duration=duration,
#             guild_snowflake=guild_snowflake,
#             is_channel_scope=is_channel_scope,
#             member_snowflake=member_snowflake,
#             message_snowflake=message_snowflake or None,
#             message_channel_snowflake=message_channel_snowflake or None,
#             reason=reason or "No reason provided.",
#             role_snowflake=role_snowflake or None,
#             target=target or None,
#         )
#
#
# async def log_xmod(
#     author_snowflake: int,
#     channel_snowflake: int,
#     display: bool,
#     guild_snowflake: int,
#     member_snowflake: int,
#     message_snowflake: int | None,
#     message_channel_snowflake: int | None,
# ):
#     duration = DurationObject(number=0, prefix="", sign=1, unit="")
#     is_channel_scope = None
#     reason = None
#     role_snowflake = None
#     target = None
#     await data_builder.save_data(
#         author_snowflake=author_snowflake or None,
#         channel_snowflake=channel_snowflake,
#         duration=duration,
#         guild_snowflake=guild_snowflake,
#         identifier="xmod",
#         member_snowflake=member_snowflake,
#         reason=reason or "No reason provided.",
#         role_snowflake=role_snowflake or None,
#         target=target or None,
#     )
#     if display:
#         await stream_service.send_log(
#             author_snowflake=author_snowflake or None,
#             channel_snowflake=channel_snowflake,
#             identifier="xmod",
#             duration=duration,
#             guild_snowflake=guild_snowflake,
#             is_channel_scope=is_channel_scope,
#             member_snowflake=member_snowflake,
#             message_snowflake=message_snowflake or None,
#             message_channel_snowflake=message_channel_snowflake or None,
#             reason=reason or "No reason provided.",
#             role_snowflake=role_snowflake or None,
#             target=target or None,
#         )
