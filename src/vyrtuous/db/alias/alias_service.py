"""action.py The purpose of this program is to be a child of DatabaseFactory and the parent to all moderation actions.

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

from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import discord

from vyrtuous.base.service import Service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.discord_object_service import DiscordObject
from vyrtuous.commands.fields.duration import DurationError, DurationObject
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.alias.alias import Alias
from vyrtuous.db.mgmt.cap.cap import Cap
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_channels,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_set_pages,
)
from vyrtuous.utils.dir_to_classes import dir_to_classes
from vyrtuous.utils.emojis import get_random_emoji


class AliasService(Service):

    lines, pages = [], []
    model = None
    ALIAS_MAP = {}

    @classmethod
    async def execute(cls, information: dict, message, state):
        kwargs = information["updated_kwargs"]
        obj = await cls.model.select(**kwargs, singular=True)
        if obj:
            if hasattr(kwargs, "role_snowflake"):
                del kwargs["role_snowflake"]
            executor_role = await PermissionService.resolve_highest_role(**kwargs)
            cap = await Cap.select(
                **kwargs, category=cls.model.identifier, singular=True
            )
            duration_seconds = DurationObject.from_expires_in(
                obj.expires_in
            ).to_seconds()
            passed = True
            if cap:
                if (
                    duration_seconds > cap.duration_seconds
                    and executor_role == "Moderator"
                ):
                    passed = False
            else:
                if (
                    duration_seconds > DurationObject("8h").to_seconds
                    and executor_role == "Moderator"
                ):
                    passed = False
            if not passed:
                raise DurationError()
            await cls.undo(information=information, message=message, state=state)
        else:
            await cls.enforce(information=information, message=message, state=state)

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        bot = DiscordBot.get_instance()
        dictionary = {}
        aliases = await Alias.select(singular=False, **where_kwargs)
        for alias in aliases:
            dictionary.setdefault(alias.guild_snowflake, {"channels": {}})
            dictionary[alias.guild_snowflake]["channels"].setdefault(
                alias.channel_snowflake, {"aliases": {}}
            )
            dictionary[alias.guild_snowflake]["channels"][alias.channel_snowflake][
                "aliases"
            ].setdefault(alias.category, {})[alias.alias_name] = []
            if alias.category == "role":
                guild = bot.get_guild(alias.guild_snowflake)
                if guild:
                    role = guild.get_role(alias.role_snowflake)
                    dictionary[alias.guild_snowflake]["channels"][
                        alias.channel_snowflake
                    ]["aliases"][alias.category][alias.alias_name] = role.mention
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_guilds:
                AliasService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_channels:
                AliasService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Command Aliases"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await AliasService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_dictionary in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                for category, alias_data in channel_dictionary["aliases"].items():
                    AliasService.lines.append(f"{category}")
                    for name, role_mention in alias_data.items():
                        if category == "role":
                            AliasService.lines.append(f"  ↳ `{name}` -> {role_mention}")
                        else:
                            AliasService.lines.append(f"  ↳ `{name}`")
                    field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(AliasService.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, AliasService.pages, title, guild.name)
                    AliasService.lines = []
                    field_count = 0
            if AliasService.lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(AliasService.lines),
                    inline=False,
                )
            AliasService.pages.append(embed)
        return AliasService.pages

    @classmethod
    async def delete_alias(cls, alias_name, default_kwargs):
        bot = DiscordBot.get_instance()
        guild_snowflake = default_kwargs.get("guild_snowflake", None)
        where_kwargs = {
            "alias_name": alias_name,
            "guild_snowflake": int(guild_snowflake),
        }
        alias = await Alias.select(singular=True, **where_kwargs)
        if not alias:
            return f"No aliases found for `{alias_name}`."
        guild = bot.get_guild(guild_snowflake)
        channel = guild.get_channel(alias.channel_snowflake)
        if getattr(alias, "role_snowflake"):
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel.mention} "
                f" and role {alias.role_mention} deleted successfully."
            )
        else:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel.mention} "
                f"deleted successfully."
            )
        await Alias.delete(**where_kwargs)
        return msg

    @classmethod
    async def create_alias(cls, alias_name, category, channel_dict, role_dict):
        where_kwargs = {}
        where_kwargs.update(channel_dict.get("columns", None))
        msg = (
            f"Alias `{alias_name}` of type `{category}` "
            f"created successfully for channel {channel_dict.get("mention", None)}."
        )
        alias = await Alias.select(category=category, **where_kwargs, singular=True)
        if alias and alias.category != "role":
            return (
                f"Alias of type `{category}` "
                f"already exists for this channel {channel_dict.get("mention", None)}."
            )
        if role_dict:
            where_kwargs.update(role_dict.get("columns", None))
            msg = (
                f"Alias `{alias_name}` of type `{category}` "
                f"created successfully for channel {channel_dict.get("mention", None)} with role {role_dict.get("mention", None)}."
            )
        alias = Alias(alias_name=alias_name, category=str(category), **where_kwargs)
        await alias.create()
        return msg

    @classmethod
    def fill_map(cls, alias, args) -> Dict[str, Tuple[int, str]]:
        map = alias.ARGS_MAP
        sorted_args = sorted(map.items(), key=lambda x: x[1])
        kwargs = {}
        for i, (key, pos) in enumerate(sorted_args):
            if i == len(sorted_args) - 1:
                value = (
                    " ".join(str(a) for a in args[pos - 1 :])
                    if len(args) >= pos
                    else ""
                )
            else:
                value = str(args[pos - 1]) if len(args) >= pos else ""
            kwargs[key] = (pos, value)
        return kwargs

    @classmethod
    async def build(cls, message):
        information = {}
        bot = DiscordBot.get_instance()
        do = DiscordObject(message=message)
        args = (
            message.content[len(bot.config["discord_command_prefix"]) :].strip().split()
        )
        alias_name = args[0]
        alias_entry = await Alias.select(
            alias_name=alias_name,
            guild_snowflake=message.guild.id,
            singular=True,
        )
        if not alias_entry:
            return
        alias_category = str(alias_entry.category)
        channel_snowflake = int(alias_entry.channel_snowflake)
        guild_snowflake = int(alias_entry.guild_snowflake)
        member_snowflake = int(message.author.id)
        default_kwargs = {
            "channel_snowflake": channel_snowflake,
            "guild_snowflake": guild_snowflake,
            "member_snowflake": member_snowflake,
        }
        information["updated_kwargs"] = default_kwargs.copy()
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous/db"))
        typed_aliases = dir_to_classes(dir_paths=dir_paths, parent=Alias)
        for alias in typed_aliases:
            if alias.category == alias_category:
                information["alias"] = alias
                break
        else:
            return
        kwargs = information["alias"].service.fill_map(
            alias=information["alias"], args=args
        )
        information["executor_role"] = await PermissionService.resolve_highest_role(
            **information["updated_kwargs"]
        )
        for field, tuple in kwargs.items():
            value = tuple[1]
            if field == "duration":
                if not value:
                    duration = DurationObject("8h")
                else:
                    duration = DurationObject(value)
                information["duration"] = duration
                cap = await Cap.select(
                    category=information["alias"].category,
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                    singular=True,
                )
                if not hasattr(cap, "duration_seconds"):
                    information["cap_duration"] = DurationObject("8h").to_seconds()
                else:
                    information["cap_duration"] = cap.duration_seconds
                if (
                    int(duration.to_timedelta().total_seconds())
                    > information["cap_duration"]
                    or duration.number == 0
                ):
                    if information["executor_role"] == "Moderator":
                        raise DurationError(
                            cap_duration=DurationObject.from_seconds(
                                information["cap_duration"]
                            ).duration
                        )
                information["expires_in"] = (
                    None
                    if duration.number == 0
                    else datetime.now(timezone.utc) + duration.to_timedelta()
                )
            elif field == "member":
                member_dict = await do.determine_from_target(target=value)
                await PermissionService.has_equal_or_lower_role(
                    updated_kwargs=default_kwargs,
                    member_snowflake=member_dict.get("id", None),
                )
                information["updated_kwargs"].update(
                    {"member_snowflake": member_dict.get("id", None)}
                )
            elif field == "reason":
                if not value:
                    information["reason"] = "No reason provided."
                else:
                    information["reason"] = value
        if getattr(alias, "role_snowflake"):
            information["updated_kwargs"].update(
                {"role_snowflake": int(alias.role_snowflake)}
            )
        return information
