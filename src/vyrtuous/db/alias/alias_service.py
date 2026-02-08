"""!/bin/python3
alias_service.py The purpose of this program is to extend Service to service aliases.

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

from vyrtuous.base.record_service import RecordService
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


class AliasService:
    lines, pages = [], []
    model = Alias

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

        alias_n = 0
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
                if not channel:
                    continue
                for category, alias_data in channel_dictionary["aliases"].items():
                    AliasService.lines.append(f"{category}")
                    for name, role_mention in alias_data.items():
                        if category == "role":
                            AliasService.lines.append(f"  ↳ `{name}` -> {role_mention}")
                        else:
                            AliasService.lines.append(f"  ↳ `{name}`")
                        alias_n += 1
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
        if AliasService.pages:
            AliasService.pages[0].description = f"**({alias_n})**"
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
            role = guild.get_role(alias.role_snowflake)
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel.mention} "
                f" and role {role.mention} deleted successfully."
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
            f"created successfully for channel {channel_dict.get('mention', None)}."
        )
        alias = await Alias.select(category=category, **where_kwargs, singular=True)
        if alias and alias.category != "role":
            return (
                f"Alias of type `{category}` "
                f"already exists for this channel {channel_dict.get('mention', None)}."
            )
        if role_dict:
            where_kwargs.update(role_dict.get("columns", None))
            msg = (
                f"Alias `{alias_name}` of type `{category}` "
                f"created successfully for channel {channel_dict.get('mention', None)} with role {role_dict.get('mention', None)}."
            )
        alias = Alias(alias_name=alias_name, category=str(category), **where_kwargs)
        await alias.create()
        return msg

    @classmethod
    def alias_category_to_alias(cls, alias_category):
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous/db"))
        typed_aliases = dir_to_classes(dir_paths=dir_paths, parent=Alias)
        for a in typed_aliases:
            if a.category == alias_category:
                return a()
        raise NotAlias()

    # await PermissionService.has_equal_or_lower_role(
    #     updated_kwargs=default_kwargs,
    #     member_snowflake=member_dict.get("id", None),
    # )
