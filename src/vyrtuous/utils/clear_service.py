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

import importlib.util
import inspect
from pathlib import Path

import discord


class ClearService:
    def __init__(
        self, *, database_factory=None, moderator_service=None, sysadmin_service=None
    ):
        self.__database_factory = database_factory
        self.__moderator_service = moderator_service
        self.__sysadmin_service = sysadmin_service

    def dir_to_classes(self, dir_paths, *, attr=None):
        classes = []
        for dir_path in dir_paths:
            for py_file in dir_path.rglob("*.py"):
                if py_file.name == "__init__.py":
                    continue
                module_name = py_file.stem
                spec = importlib.util.spec_from_file_location(module_name, str(py_file))
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                for _, cls in inspect.getmembers(module, inspect.isclass):
                    if cls.__module__ != module.__name__:
                        continue
                    if getattr(cls, "__skip_db_discovery__", False):
                        continue
                    if attr in getattr(cls, "__annotations__", {}) or not attr:
                        classes.append(cls)
        return classes

    async def clear(
        self,
        category,
        context,
        object_dict,
        target,
        view,
        where_kwargs,
    ):
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous"))
        if isinstance(object_dict.get("object", None), discord.Member):
            await self.__moderator_service.has_equal_or_lower_role(
                **context.to_dict(),
                target_member_snowflake=object_dict.get("id", None),
            )
            if view.result:
                for obj in self.dir_to_classes(
                    dir_paths=dir_paths, attr="member_snowflake"
                ):
                    if str(category) == "all":
                        await self.__database_factory.delete_by_cls(obj, **where_kwargs)
                        msg = f"Deleted all associated database information for {object_dict.get('mention', None)}."
                    elif str(category).lower() == obj.identifier:
                        await self.__database_factory.delete_by_cls(obj, **where_kwargs)
                        msg = f"Deleted all associated {category} records for {object_dict.get('mention', None)}."
        elif isinstance(object_dict.get("object", None), discord.abc.GuildChannel):
            await self.__moderator_service.check(
                **object_dict.get("columns", None), lowest_role="Guild Owner"
            )
            if view.result:
                for obj in self.dir_to_classes(
                    dir_paths=dir_paths, attr="channel_snowflake"
                ):
                    if category == "all":
                        await self.__database_factory.delete_by_cls(obj, **where_kwargs)
                        msg = f"Deleted all database information for {object_dict.get('mention')}."
                    elif str(category).lower() == obj.identifier:
                        await self.__database_factory.delete_by_cls(obj, **where_kwargs)
                        msg = f"Deleted all associated {category} records in {object_dict.get('mention', None)}."
        elif isinstance(object_dict.get("object", None), discord.Guild):
            await self.__moderator_service.check(
                **object_dict.get("columns", None), lowest_role="Developer"
            )
            if view.result:
                attributes = [
                    "channel_snowflake",
                    "guild_snowflake",
                    "member_snowflake",
                ]
                objects = [
                    obj
                    for attr in attributes
                    for obj in self.dir_to_classes(dir_paths=dir_paths, attr=attr)
                ]
                for obj in objects:
                    if str(category).lower() == "all":
                        await self.__database_factory.delete_by_cls(obj, **where_kwargs)
                        msg = f"Deleted all database information for {object_dict.get('name')}."
                    elif str(category).lower() == obj.identifier:
                        await self.__database_factory.delete_by_cls(obj, **where_kwargs)
                        msg = f"Deleted all associated {category} in {object_dict.get('name', None)}."
        elif target == "all" and await self.__sysadmin_service.is_sysadmin(
            member_snowflake=context.to_dict().get("member_snowflake")
        ):
            if view.result:
                for obj in self.dir_to_classes(dir_paths=dir_paths, attr=None):
                    await self.__database_factory.delete_by_cls(obj, **where_kwargs)
                    msg = "Deleted all database entries."
        else:
            msg = f"Invalid target ({target})."
        return msg
