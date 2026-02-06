"""!/bin/python3
clear_service.py The purpose of this program is to service the clear command.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from pathlib import Path

import discord

from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.base.service import Service
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.roles.admin.administrator import AdministratorRole
from vyrtuous.db.roles.sysadmin.sysadmin_service import is_sysadmin
from vyrtuous.utils.dir_to_classes import dir_to_classes


class ClearService(Service):

    @classmethod
    async def clear(
        cls,
        category,
        default_kwargs,
        object_dict,
        target,
        view,
        where_kwargs,
    ):
        guild_snowflake = default_kwargs.get("guild_snowflake", None)
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous/db"))
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(object_dict.get("columns", None))
        if isinstance(object_dict.get("object", None), discord.Member):
            await PermissionService.has_equal_or_lower_role(
                updated_kwargs=default_kwargs,
                member_snowflake=object_dict.get("id", None),
            )
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths, parent=DatabaseFactory):
                    if "member_snowflake" in getattr(obj, "__annotations__", {}):
                        if str(category) == "all":
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated database information for {object_dict.get('mention', None)}."
                        elif str(category).lower() == obj.identifier:
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated {category} records for {object_dict.get('mention', None)}."
        elif isinstance(object_dict.get("object", None), discord.abc.GuildChannel):
            await PermissionService.check(
                updated_kwargs=updated_kwargs, lowest_role="Guild Owner"
            )
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths, parent=DatabaseFactory):
                    if "channel_snowflake" in getattr(obj, "__annotations__", {}):
                        if category == "all":
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all database information for {object_dict.get('mention')}."
                        elif str(category).lower() == obj.identifier:
                            await obj.delete(**where_kwargs)
                            msg = f"Deleted all associated {category} records in {object_dict.get('mention', None)}."
        elif isinstance(object_dict.get("object", None), discord.Guild):
            await PermissionService.check(
                updated_kwargs=updated_kwargs, lowest_role="Developer"
            )
            if view.result:
                for obj in dir_to_classes(dir_paths=dir_paths, parent=DatabaseFactory):
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
                        elif str(category).lower() == obj.identifier:
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
            member_snowflake=default_kwargs.get("member_snowflake")
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
