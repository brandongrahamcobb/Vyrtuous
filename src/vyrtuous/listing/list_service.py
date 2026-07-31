"""!/bin/python3
list_service.py The purpose of this program is to service listing.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Union

import discord

from vyrtuous.bot.discord_bot import DiscordBot

CHUNK_SIZE = 12

from typing import TypeVar


@dataclass(frozen=True)
class AutoAssignRoleDictionary:
    data: Dict[int, Dict[str, Dict[int, dict]]] = field(default_factory=dict)


@dataclass(frozen=True)
class AliasDictionary:
    data: Dict[
        int, Dict[str, Dict[int, Dict[str, Dict[str, Dict[str, list | str]]]]]
    ] = field(default_factory=dict)


@dataclass
class AutoMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, Any]]]]] = field(
        default_factory=dict
    )


@dataclass
class BanDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )


@dataclass
class CapDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, int]]]]] = field(
        default_factory=dict
    )


@dataclass
class FlagDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, str]]]]]] = field(
        default_factory=dict
    )


@dataclass(frozen=True)
class HeroDictionary:
    data: Dict[int, Dict[str, Dict[int, bool]]] = field(default_factory=dict)


@dataclass(frozen=True)
class PermissionDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, List[str]]]]] = field(
        default_factory=dict
    )


@dataclass
class ServerMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, dict]]]] = field(default_factory=dict)


@dataclass
class StreamDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, List[int]]]]] = field(
        default_factory=dict
    )


@dataclass
class TextMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )


@dataclass
class VeganDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, dict]]]] = field(default_factory=dict)


@dataclass
class VideoChannelDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict]]] = field(default_factory=dict)


@dataclass
class VoiceMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )


DictT = TypeVar(
    "DictT",
    bound=Union[
        BanDictionary,
        AutoAssignRoleDictionary,
        AliasDictionary,
        AutoMuteDictionary,
        CapDictionary,
        FlagDictionary,
        HeroDictionary,
        PermissionDictionary,
        ServerMuteDictionary,
        StreamDictionary,
        TextMuteDictionary,
        VeganDictionary,
        VideoChannelDictionary,
        VoiceMuteDictionary,
    ],
)


async def process_dictionary(cls, dictionary) -> DictT:
    skipped_pages = {}
    data = clean_dictionary(
        dictionary=dictionary,
    )
    kwargs: dict[str, object] = {"data": data}
    return cls(**kwargs)


def generate_skipped_set_pages(skipped, title) -> list[discord.Embed]:
    field_count = 0
    pages: list[discord.Embed] = []
    embed = discord.Embed(title=title, description="\u200b", color=discord.Color.blue())
    lines: list[str] = []
    for snowflake in skipped:
        if field_count >= CHUNK_SIZE:
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


def generate_skipped_dict_pages(skipped, title) -> list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    field_count = 0
    pages: list[discord.Embed] = []
    for guild_snowflake, data in skipped.items():
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            color=discord.Color.red(), title=f"{title} ({guild.name})"
        )
        field_count = 0
        lines: list[str] = []
        for snowflake in data:
            if field_count >= CHUNK_SIZE:
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


def generate_skipped_guilds(dictionary: dict) -> set[int]:
    bot: DiscordBot = DiscordBot.get_instance()
    skipped_guilds = set()
    for guild_snowflake in dictionary:
        if not bot.get_guild(guild_snowflake):
            skipped_guilds.add(guild_snowflake)
    return skipped_guilds


def generate_skipped_channels(dictionary: dict) -> dict[int, list[int]]:
    bot: DiscordBot = DiscordBot.get_instance()
    skipped_channels: dict[int, list[int]] = {}
    for guild_snowflake, guild_data in dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        for channel_snowflake in guild_data.get("channels", {}):
            if not guild.get_channel(channel_snowflake):
                skipped_channels.setdefault(guild_snowflake, []).append(
                    channel_snowflake
                )
    return skipped_channels


def generate_skipped_members(dictionary: dict) -> dict[int, list[int]]:
    bot: DiscordBot = DiscordBot.get_instance()
    skipped_members: dict[int, list[int]] = {}
    for guild_snowflake, guild_data in dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        for member_snowflake in guild_data.get("members", {}):
            if not guild.get_member(member_snowflake):
                skipped_members.setdefault(guild_snowflake, []).append(member_snowflake)
    return skipped_members


def generate_skipped_roles(dictionary: dict) -> dict[int, list[int]]:
    bot: DiscordBot = DiscordBot.get_instance()
    skipped_roles: dict[int, list[int]] = {}
    for guild_snowflake, guild_data in dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        for role_snowflake in guild_data.get("roles", {}):
            if not guild.get_role(role_snowflake):
                skipped_roles.setdefault(guild_snowflake, []).append(role_snowflake)
    return skipped_roles


def clean_dictionary(
    dictionary: dict,
) -> dict[Any, dict[str, dict[Any, Any]]]:
    cleaned = {}
    for guild_snowflake, guild_data in dictionary.items():
        channels = {
            channel_snowflake: entries
            for channel_snowflake, entries in guild_data.get("channels", {}).items()
        }
        members = {
            member_snowflake: entries
            for member_snowflake, entries in guild_data.get("members", {}).items()
        }
        messages = {
            message_snowflake: entries
            for message_snowflake, entries in guild_data.get("messages", {}).items()
        }
        roles = {
            role_snowflake: entries
            for role_snowflake, entries in guild_data.get("roles", {}).items()
        }
        snowflakes = {
            snowflake: entries
            for snowflake, entries in guild_data.get("snowflakes", {}).items()
        }
        cleaned[guild_snowflake] = {
            "channels": channels,
            "members": members,
            "messages": messages,
            "roles": roles,
            "snowflakes": snowflakes,
        }
    return cleaned


async def generate_skipped_messages(dictionary: dict) -> dict[int, list[int]]:
    bot: DiscordBot = DiscordBot.get_instance()
    skipped_messages: dict[int, list[int]] = {}
    for guild_snowflake, guild_data in dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        for channel_snowflake, channel_logs in guild_data.get("channels", {}).items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            if not isinstance(channel, (discord.VoiceChannel, discord.TextChannel)):
                continue
            members = channel_logs.get("members")
            if not isinstance(members, dict):
                continue
            for _, member_data in members.items():
                if not isinstance(member_data, dict):
                    continue
                messages = member_data.get("messages")
                if not isinstance(messages, dict):
                    continue
                for message_snowflake, _ in messages.items():
                    try:
                        await channel.fetch_message(message_snowflake)
                    except Exception:
                        skipped_messages.setdefault(guild_snowflake, []).append(
                            message_snowflake
                        )
    return skipped_messages


def flush_page(embed, pages, title, guild_name) -> discord.Embed:
    pages.append(embed)
    return discord.Embed(
        title=title,
        description=f"{guild_name} continued...",
        color=discord.Color.blue(),
    )
