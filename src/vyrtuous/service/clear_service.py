from pathlib import Path

import discord

from vyrtuous.db.roles.administrator import AdministratorRole
from vyrtuous.service.roles.sysadmin_service import is_sysadmin
from vyrtuous.utils.check import check, has_equal_or_lower_role
from vyrtuous.utils.dir_to_classes import dir_to_classes


class ClearService:

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
            await has_equal_or_lower_role(
                snowflake_kwargs=snowflake_kwargs,
                member_snowflake=object_dict.get("id", None),
            )
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if "member_snowflake" in getattr(obj, "__annotations__", {}):
                        if str(category) == "all":
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated database information for {object_dict.get('mention', None)}."
                        elif str(category).lower() == obj.category:
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated {category} records for {object_dict.get('mention', None)}."
        elif isinstance(object_dict.get("object", None), discord.abc.GuildChannel):
            await check(snowflake_kwargs=snowflake_kwargs, lowest_role="Guild Owner")
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if "channel_snowflake" in getattr(obj, "__annotations__", {}):
                        if category == "all":
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all database information for {object_dict.get('mention')}."
                        elif str(category).lower() == obj.category:
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated {category} records in {object_dict.get('mention', None)}."
        elif isinstance(object_dict.get("object", None), discord.Guild):
            await check(snowflake_kwargs=snowflake_kwargs, lowest_role="Developer")
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths):
                    if any(
                        getattr(obj, attr, None) is not None
                        or attr in getattr(obj, "__annotations__", {})
                        for attr in (
                            "channel_snowflake",
                            "guild_snowflake",
                            "member_snowflake",
                        )
                    ):
                        if str(category).lower() == "all":
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all database information for {object_dict.get('name')}."
                        elif str(category).lower() == obj.category:
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated {category} in {object_dict.get('name', None)}."
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
