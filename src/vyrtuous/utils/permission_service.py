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

import discord


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

    async def build_clean_dictionary(self, channel_objs, is_at_home, me):
        dictionary = {}
        pages = []
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
        skipped_channels = self.__dictionary_service.generate_skipped_channels(
            dictionary
        )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_channels:
                pages.extend(
                    self.__dictionary_service.generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
            if skipped_guilds:
                pages.extend(
                    self.__dictionary_service.generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    async def build_pages(self, is_at_home, channel_objs, default_kwargs):
        lines, pages = []
        channel_snowflake = default_kwargs.get("channel_snowflake", None)
        guild_snowflake = default_kwargs.get("guild_snowflake", None)
        guild = self.__bot.get_guild(guild_snowflake)
        title = f"{self.__emoji.get_random_emoji()} {self.__bot.user.display_name} Missing Permissions"

        dictionary = await self.build_clean_dictionary(
            channel_objs=channel_objs, is_at_home=is_at_home, me=guild.me
        )

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
