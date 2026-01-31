from pathlib import Path

import discord

from vyrtuous.utils.dir_to_classes import dir_to_classes
from vyrtuous.db.roles.administrator import (
    AdministratorRole,
)
from vyrtuous.db.roles.sysadmin import is_sysadmin
from vyrtuous.utils.check import (
    check,
    has_equal_or_lower_role,
)


class Clear:

    @classmethod
    async def clear(
        cls,
        category,
        where_kwargs,
        object_dict,
        snowflake_kwargs,
        target,
        view,
    ):
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        dir_paths = []
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/infractions")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/mgmt")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/roles")
        dir_paths.append(Path(__file__).resolve().parents[1] / "db/rooms")
        if isinstance(object_dict.get("object", None), discord.Member):
            has_equal_or_lower_role(
                snowflake_kwargs=snowflake_kwargs,
                member_snowflake=object_dict.get("id", None),
            )
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if "member" in obj.SCOPES:
                        if str(category) == "all":
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated database information for {object_dict.get('mention', None)}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} for {object_dict.get('mention', None)}."
        elif isinstance(object_dict.get("object", None), discord.abc.GuildChannel):
            await check(snowflake_kwargs=snowflake_kwargs, lowest_role="Guild Owner")
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if "channel" in obj.SCOPES:
                        if category == "all":
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all database information for {object_dict.get('mention')}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} in {object_dict.get('mention', None)}."
        elif isinstance(object_dict.get("object", None), discord.Guild):
            await check(snowflake_kwargs=snowflake_kwargs, lowest_role="Guild Owner")
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if any(
                        scope in obj.SCOPES for scope in ("guild", "channel", "member")
                    ):
                        if str(category).lower() == "all":
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all database information for {object_dict.get('name')}."
                        elif str(category).lower() == obj.CATEGORY:
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated {obj.PLURAL.lower()} in {object_dict.get('name', None)}."
                        elif isinstance(obj, AdministratorRole):
                            administrator_roles = AdministratorRole.select(
                                guild_snowflake=int(guild_snowflake),
                            )
                            for administrator_role in administrator_roles:
                                await obj.revoke_role(
                                    guild_snowflake=int(guild_snowflake),
                                    role_snowflake=administrator_role.role_snowflake,
                                )
        elif target == "all" and await is_sysadmin(
            member_snowflake=snowflake_kwargs.get("member_snowflake")
        ):
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    await obj.delete(**where_kwargs)
                    msg = "Deleted all database entries."
                    if isinstance(obj, AdministratorRole):
                        administrator_roles = AdministratorRole.select(
                            guild_snowflake=int(guild_snowflake),
                        )
                        for administrator_role in administrator_roles:
                            await obj.revoke_role(
                                guild_snowflake=int(guild_snowflake),
                                role_snowflake=administrator_role.role_snowflake,
                            )
        else:
            msg = f"Invalid target ({target})."
        return msg
