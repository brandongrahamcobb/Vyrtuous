"""!/bin/python3
permission_service.py The purpose of this program is to provide the service for deciding whether a member has sufficient permissions.

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

from dataclasses import dataclass, field
from typing import Dict, List

import discord


@dataclass(frozen=True)
class PermissionDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, List[str]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guild: List[discord.Embed] = field(default_factory=list)


class PermissionService:
    __CHUNK_SIZE = 7
    __TARGET_PERMISSIONS = (
        "add_reactions",
        "manage_messages",
        "move_members",
        "mute_members",
        "send_messages",
        "view_channel",
    )

    def __init__(self, *, bot=None, dictionary_service=None, emoji=None):
        self.__bot = bot
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji

    def build_dictionary(self, channel_objs, me):
        dictionary = {}
        for channel in channel_objs:
            permissions = channel.permissions_for(me)
            missing = []
            for permission in self.__TARGET_PERMISSIONS:
                if not getattr(permissions, permission):
                    missing.append(permission)
            if not missing:
                continue
            dictionary.setdefault(channel.guild.id, {"channels": {}})
            dictionary[channel.guild.id]["channels"].setdefault(channel.id, {})
            dictionary[channel.guild.id]["channels"][channel.id].update(
                {"permissions": missing}
            )
        return dictionary

    async def build_pages(self, channel_objs, context, is_at_home):
        lines, pages = []
        guild = self.__bot.get_guild(context.guild.id)
        title = f"{self.__emoji.get_random_emoji()} {self.__bot.user.display_name} Missing Permissions"

        dictionary = self.build_dictionary(channel_objs=channel_objs, me=guild.me)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=PermissionDictionary, dictionary=dictionary
        )
        if is_at_home:
            pages.extend(processed_dictionary.skipped_channels)
            pages.extend(processed_dictionary.skipped_guilds)

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(f"Channel: {channel.mention}")
                for section_name, permissions in channel_data.items():
                    for permission in permissions:
                        lines.append(f"  ↳ {permission}")
                field_count += 1
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                    lines = []
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)
        return pages
