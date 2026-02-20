"""!/bin/python3

dictionary.py The purpose of this program is to manage list command logic.

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


class DictionaryService:
    __CHUNK_SIZE = 7

    def __init__(self, *, bot=None):
        self.__bot = bot

    async def process_dictionary(self, cls, dictionary):
        skipped_pages = {}
        skipped_channels = self.generate_skipped_channels(dictionary)
        skipped_guilds = self.generate_skipped_guilds(dictionary)
        skipped_members = self.generate_skipped_members(dictionary)
        skipped_messages = await self.generate_skipped_messages(dictionary)
        skipped_roles = self.generate_skipped_roles(dictionary)
        data = self.clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
            skipped_messages=skipped_messages,
            skipped_roles=skipped_roles,
        )
        if hasattr(cls, "skipped_channels") and cls.skipped_channels:
            skipped_pages["skipped_channels"] = self.generate_skipped_dict_pages(
                skipped=cls.skipped_channels,
                title="Skipped Channels in Server",
            )
        if hasattr(cls, "skipped_guilds") and cls.skipped_guilds:
            skipped_pages["skipped_guilds"] = self.generate_skipped_set_pages(
                skipped=cls.skipped_guilds,
                title="Skipped Servers",
            )
        if hasattr(cls, "skipped_members") and cls.skipped_members:
            skipped_pages["skipped_members"] = self.generate_skipped_dict_pages(
                skipped=cls.skipped_members,
                title="Skipped Members in Server",
            )
        if hasattr(cls, "skipped_messages") and cls.skipped_messages:
            skipped_pages["skipped_messages"] = self.generate_skipped_dict_pages(
                skipped=cls.skipped_messages,
                title="Skipped Messages in Server",
            )
        if hasattr(cls, "skipped_roles") and cls.skipped_roles:
            skipped_pages["skipped_roles"] = self.generate_skipped_dict_pages(
                skipped=cls.skipped_roles,
                title="Skipped Roles in Server",
            )
        return cls(data=data, **skipped_pages)

    def generate_skipped_set_pages(self, skipped, title):
        field_count = 0
        pages = []
        embed = discord.Embed(
            title=title, description="\u200b", color=discord.Color.blue()
        )
        lines = []
        for snowflake in skipped:
            if field_count >= self.__CHUNK_SIZE:
                embed.description = "\n".join(lines)
                pages.append(embed)
                embed = discord.Embed(
                    title=f"{title} continued...", color=discord.Color.red()
                )
                lines = []
                field_count = 0
            lines.append(str(snowflake))
            field_count += 1
        embed.description = "\n".join(lines)
        pages.append(embed)
        return pages

    def generate_skipped_dict_pages(self, skipped, title):
        field_count = 0
        pages = []
        for guild_snowflake, list in skipped.items():
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                color=discord.Color.red(), title=f"{title} ({guild.name})"
            )
            field_count = 0
            lines = []
            for snowflake in list:
                if field_count >= self.__CHUNK_SIZE:
                    embed.description = "\n".join(lines)
                    pages.append(embed)
                    embed = discord.Embed(
                        color=discord.Color.red(),
                        title=f"{title} ({guild_snowflake}) continued...",
                    )
                    field_count = 0
                    lines = []
                lines.append(str(snowflake))
                field_count += 1
            embed.description = "\n".join(lines)
            pages.append(embed)
        return pages

    def generate_skipped_guilds(self, dictionary: dict) -> set:
        skipped_guilds = set()
        for guild_snowflake in dictionary:
            if not self.__bot.get_guild(guild_snowflake):
                skipped_guilds.add(guild_snowflake)
        return skipped_guilds

    def generate_skipped_channels(self, dictionary: dict) -> dict:
        skipped_channels = {}
        for guild_snowflake, guild_data in dictionary.items():
            guild = self.__bot.get_guild(guild_snowflake)
            if not guild:
                continue
            for channel_snowflake in guild_data.get("channels", {}):
                if not guild.get_channel(channel_snowflake):
                    skipped_channels.setdefault(guild_snowflake, []).append(
                        channel_snowflake
                    )
        return skipped_channels

    def generate_skipped_members(self, dictionary: dict) -> dict:
        skipped_members = {}
        for guild_snowflake, guild_data in dictionary.items():
            guild = self.__bot.get_guild(guild_snowflake)
            if not guild:
                continue
            for member_snowflake in guild_data.get("members", {}):
                if not guild.get_member(member_snowflake):
                    skipped_members.setdefault(guild_snowflake, []).append(
                        member_snowflake
                    )
        return skipped_members

    def generate_skipped_roles(self, dictionary: dict) -> dict:
        skipped_roles = {}
        for guild_snowflake, guild_data in dictionary.items():
            guild = self.__bot.get_guild(guild_snowflake)
            if not guild:
                continue
            for role_snowflake in guild_data.get("roles", {}):
                if not guild.get_role(role_snowflake):
                    skipped_roles.setdefault(guild_snowflake, []).append(role_snowflake)
        return skipped_roles

    def clean_dictionary(
        self,
        dictionary: dict,
        *,
        skipped_guilds: set | None = None,
        skipped_channels: dict | None = None,
        skipped_members: dict | None = None,
        skipped_messages: dict | None = None,
        skipped_roles: dict | None = None,
        skipped_snowflakes: dict | None = None,
    ) -> dict:
        cleaned = {}
        skipped_guilds = skipped_guilds or set()
        skipped_channels = skipped_channels or {}
        skipped_members = skipped_members or {}
        skipped_messages = skipped_messages or {}
        skipped_roles = skipped_roles or {}
        skipped_snowflakes = skipped_snowflakes or {}
        for guild_snowflake, guild_data in dictionary.items():
            if guild_snowflake in skipped_guilds:
                continue
            channels = {
                channel_snowflake: entries
                for channel_snowflake, entries in guild_data.get("channels", {}).items()
                if channel_snowflake not in skipped_channels.get(guild_snowflake, [])
            }
            members = {
                member_snowflake: entries
                for member_snowflake, entries in guild_data.get("members", {}).items()
                if member_snowflake not in skipped_members.get(guild_snowflake, [])
            }
            messages = {
                message_snowflake: entries
                for message_snowflake, entries in guild_data.get("messages", {}).items()
                if message_snowflake not in skipped_messages.get(guild_snowflake, [])
            }
            roles = {
                role_snowflake: entries
                for role_snowflake, entries in guild_data.get("roles", {}).items()
                if role_snowflake not in skipped_roles.get(guild_snowflake, [])
            }
            snowflakes = {
                snowflake: entries
                for snowflake, entries in guild_data.get("snowflakes", {}).items()
                if snowflake not in skipped_snowflakes.get(guild_snowflake, [])
            }
            cleaned[guild_snowflake] = {
                "channels": channels,
                "members": members,
                "messages": messages,
                "roles": roles,
                "snowflakes": snowflakes,
            }
        return cleaned

    async def generate_skipped_messages(self, dictionary: dict) -> dict:
        skipped_messages = {}
        for guild_snowflake, guild_data in dictionary.items():
            guild = self.__bot.get_guild(guild_snowflake)
            if not guild:
                continue
            for channel_snowflake, channel_logs in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                members = channel_logs.get("members")
                if not isinstance(members, dict):
                    continue
                for member_id, member_data in members.items():
                    if not isinstance(member_data, dict):
                        continue
                    messages = member_data.get("messages")
                    if not isinstance(messages, dict):
                        continue
                    for message_snowflake, message_data in messages.items():
                        try:
                            await channel.fetch_message(message_snowflake)
                        except Exception:
                            skipped_messages.setdefault(guild_snowflake, []).append(
                                message_snowflake
                            )
        return skipped_messages

    def flush_page(self, embed, pages, title, guild_name):
        pages.append(embed)
        return discord.Embed(
            title=title,
            description=f"{guild_name} continued...",
            color=discord.Color.blue(),
        )
