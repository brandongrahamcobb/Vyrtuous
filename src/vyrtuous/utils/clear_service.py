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
        self,
        *,
        ban_service=None,
        database_factory=None,
        flag_service=None,
        moderator_service=None,
        sysadmin_service=None,
        text_mute_service=None,
        voice_mute_service=None,
    ):
        self.__ban_service = ban_service
        self.__database_factory = database_factory
        self.__flag_service = flag_service
        self.__moderator_service = moderator_service
        self.__sysadmin_service = sysadmin_service
        self.__text_mute_service = text_mute_service
        self.__voice_mute_service = voice_mute_service

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
        self, category, default_ctx, object_dict, target, view, where_kwargs, source
    ):
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous"))
        if isinstance(object_dict.get("object", None), discord.Member):
            await self.__moderator_service.has_equal_or_lower_role(
                **default_ctx.to_dict(),
                target_member_snowflake=object_dict.get("id", None),
            )
            if view.result:
                for cls in self.dir_to_classes(
                    dir_paths=dir_paths, attr="member_snowflake"
                ):
                    if str(category) == "all":
                        msg = f"Deleted all associated database information for {object_dict.get('mention', None)}."
                    elif str(category).lower() == cls.identifier:
                        msg = f"Deleted all associated {category} records for {object_dict.get('mention', None)}."
                    if cls.__name__ == "Ban":
                        await self.__ban_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "Flag":
                        await self.__flag_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "TextMute":
                        await self.__text_mute_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "VoiceMute":
                        await self.__voice_mute_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    else:
                        await self.__database_factory.delete_by_cls(cls, **where_kwargs)
        elif isinstance(object_dict.get("object", None), discord.abc.GuildChannel):
            await self.__moderator_service.has_equal_or_lower_role(
                **object_dict.get("columns", None), lowest_role="Guild Owner"
            )
            if view.result:
                for cls in self.dir_to_classes(
                    dir_paths=dir_paths, attr="channel_snowflake"
                ):
                    if category == "all":
                        msg = f"Deleted all database information for {object_dict.get('mention')}."
                    elif str(category).lower() == cls.identifier:
                        msg = f"Deleted all associated {category} records in {object_dict.get('mention', None)}."
                    if cls.__name__ == "Ban":
                        await self.__ban_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "Flag":
                        await self.__flag_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "TextMute":
                        await self.__text_mute_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "VoiceMute":
                        await self.__voice_mute_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    else:
                        await self.__database_factory.delete_by_cls(cls, **where_kwargs)
        elif isinstance(object_dict.get("object", None), discord.Guild):
            await self.__moderator_service.has_equal_or_lower_role(
                **object_dict.get("columns", None), lowest_role="Developer"
            )
            if view.result:
                attributes = [
                    "channel_snowflake",
                    "guild_snowflake",
                    "member_snowflake",
                ]
                classes = [
                    cls
                    for attr in attributes
                    for cls in self.dir_to_classes(dir_paths=dir_paths, attr=attr)
                ]
                for cls in classes:
                    if str(category).lower() == "all":
                        msg = f"Deleted all database information for {object_dict.get('name')}."
                    elif str(category).lower() == obj.identifier:
                        msg = f"Deleted all associated {category} in {object_dict.get('name', None)}."
                    if cls.__name__ == "Ban":
                        await self.__ban_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "Flag":
                        await self.__flag_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "TextMute":
                        await self.__text_mute_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    elif cls.__name__ == "VoiceMute":
                        await self.__voice_mute_service.delete(
                            author=default_ctx.author,
                            kwargs=where_kwargs,
                            source=source,
                        )
                    else:
                        await self.__database_factory.delete_by_cls(cls, **where_kwargs)
        elif target == "all" and await self.__sysadmin_service.is_sysadmin(
            member_snowflake=default_ctx.to_dict().get("member_snowflake")
        ):
            if view.result:
                for cls in self.dir_to_classes(dir_paths=dir_paths, attr=None):
                    if cls.__name__ == "Ban":
                        await self.__ban_service.delete(
                            author=default_ctx.author, **where_kwargs
                        )
                    elif cls.__name__ == "Flag":
                        await self.__flag_service.delete(
                            author=default_ctx.author, **where_kwargs
                        )
                    elif cls.__name__ == "TextMute":
                        await self.__text_mute_service.delete(
                            author=default_ctx.author, **where_kwargs
                        )
                    elif cls.__name__ == "VoiceMute":
                        await self.__voice_mute_service.delete(
                            author=default_ctx.author, **where_kwargs
                        )
                    else:
                        await self.__database_factory.delete_by_cls(cls, **where_kwargs)
                    msg = "Deleted all database entries."

        else:
            msg = f"Invalid target ({target})."
        return msg
