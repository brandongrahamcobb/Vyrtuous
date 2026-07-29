from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState


def exists(permission: str) -> bool:
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    if permission.endswith(".*"):
        permission = permission[:-2]
    return permission in permission_state.nodes


def validate_groups() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    sysadmin_groups = []
    guild_owner_groups = []
    for group in permission_state.groups.values():
        if group.is_sysadmin:
            if group not in sysadmin_groups:
                sysadmin_groups.append(group)
            if len(sysadmin_groups) > 1:
                group_names = [group.name for group in sysadmin_groups]
                raise ValueError(
                    f"Multiple groups ({', '.join(group_names)}) have sysadmin privileges. Only one group can have sysadmin permissions."
                )
        if group.is_guild_owner:
            if group not in guild_owner_groups:
                guild_owner_groups.append(group)
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
            if not exists(node):
                raise ValueError(
                    f'Group "{group.name}" has unknown permission "{permission}"'
                )
        for parent in group.inheritance:
            if parent not in permission_state.groups:
                raise ValueError(
                    f'Group "{group.name}" inherits unknown group "{parent}"'
                )
