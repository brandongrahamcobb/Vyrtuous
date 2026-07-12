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

from dataclasses import dataclass, field, fields
from typing import Any, Dict, List, Union

import discord

from vyrtuous.bot.discord_bot import DiscordBot

CHUNK_SIZE = 12

from typing import TypeVar, cast


@dataclass(frozen=True)
class AdministratorRoleDictionary:
    data: Dict[int, Dict[str, Dict[int, dict]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_roles: List[discord.Embed] = field(default_factory=list)


@dataclass(frozen=True)
class AdministratorDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, bool]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass(frozen=True)
class AliasDictionary:
    data: Dict[
        int, Dict[str, Dict[int, Dict[str, Dict[str, Dict[str, list | str]]]]]
    ] = field(default_factory=dict)
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


@dataclass
class AutoMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, Any]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


@dataclass
class BanDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass
class CapDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, int]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


@dataclass
class CoordinatorDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, str]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass
class FlagDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, str]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass(frozen=True)
class HeroDictionary:
    data: Dict[int, Dict[str, Dict[int, bool]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass
class ModeratorDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, str]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass(frozen=True)
class PermissionDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, List[str]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


@dataclass
class ServerMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, dict]]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass
class StreamDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, List[int]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


@dataclass
class TextMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass
class VeganDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, dict]]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


@dataclass
class VideoChannelDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict]]] = field(default_factory=dict)
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


@dataclass
class VoiceMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


DictT = TypeVar(
    "DictT",
    bound=Union[
        BanDictionary,
        AdministratorDictionary,
        AdministratorRoleDictionary,
        AliasDictionary,
        AutoMuteDictionary,
        CapDictionary,
        CoordinatorDictionary,
        FlagDictionary,
        HeroDictionary,
        ModeratorDictionary,
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
    skipped_channels = generate_skipped_channels(dictionary)
    skipped_guilds = generate_skipped_guilds(dictionary)
    skipped_members = generate_skipped_members(dictionary)
    skipped_messages = await generate_skipped_messages(dictionary)
    skipped_roles = generate_skipped_roles(dictionary)
    data = clean_dictionary(
        dictionary=dictionary,
        skipped_channels=skipped_channels,
        skipped_guilds=skipped_guilds,
        skipped_members=skipped_members,
        skipped_messages=skipped_messages,
        skipped_roles=skipped_roles,
    )
    field_names = {f.name for f in fields(cls)}
    kwargs: dict[str, object] = {"data": data}
    if "skipped_channels" in field_names:
        kwargs["skipped_channels"] = generate_skipped_dict_pages(
            skipped=skipped_channels, title="Skipped Channels in Server"
        )
    if "skipped_guilds" in field_names:
        kwargs["skipped_guilds"] = generate_skipped_set_pages(
            skipped=skipped_guilds, title="Skipped Servers"
        )
    if "skipped_members" in field_names:
        kwargs["skipped_members"] = generate_skipped_dict_pages(
            skipped=skipped_members, title="Skipped Members in Server"
        )
    if "skipped_messages" in field_names:
        kwargs["skipped_messages"] = generate_skipped_dict_pages(
            skipped=skipped_messages, title="Skipped Messages in Server"
        )
    if "skipped_roles" in field_names:
        kwargs["skipped_roles"] = generate_skipped_dict_pages(
            skipped=skipped_roles, title="Skipped Roles in Server"
        )
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
    *,
    skipped_guilds: set | None = None,
    skipped_channels: dict | None = None,
    skipped_members: dict | None = None,
    skipped_messages: dict | None = None,
    skipped_roles: dict | None = None,
    skipped_snowflakes: dict | None = None,
) -> dict[Any, dict[str, dict[Any, Any]]]:
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
