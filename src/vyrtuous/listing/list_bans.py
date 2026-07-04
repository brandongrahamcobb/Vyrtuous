"""!/bin/python3
ban_service.py The purpose of this program is to extend AliasService to service ban infractions.

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
from typing import Any, Dict, List

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.ban import Ban
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis

MODEL = Ban


@dataclass
class BanDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


async def build_dictionary(obj) -> dict:
    database_factory = DatabaseFactory(MODEL)
    bans = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        bans = await database_factory.select(guild_snowflake=obj.id, singular=False)
    elif isinstance(obj, discord.abc.GuildChannel):
        bans = await database_factory.select(channel_snowflake=obj.id, singular=False)
    elif isinstance(obj, discord.Member):
        bans = await database_factory.select(member_snowflake=obj.id, singular=False)
    else:
        bans = await database_factory.select(singular=False)
    if bans:
        for ban in bans:
            dictionary.setdefault(ban.guild_snowflake, {"members": {}})
            dictionary[ban.guild_snowflake]["members"].setdefault(
                ban.member_snowflake, {"bans": {}}
            )
            dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                "bans"
            ].setdefault(ban.channel_snowflake, {})
            dictionary[ban.guild_snowflake]["members"][ban.member_snowflake]["bans"][
                ban.channel_snowflake
            ] = {
                "reason": ban.reason,
                "expires_in": ban.expires_in,
                "blacklisted": ban.blacklisted,
            }
    return dictionary


async def build_pages(is_at_home: bool, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    duration_builder = DurationBuilder()
    lines: list[str] = []
    pages: list[discord.Embed] = []
    thumbnail = False

    obj_name = "All Servers"
    if not isinstance(obj, int):
        obj_name = obj.name
    else:
        simplified_member = bot.registry.get(MemberState).active.get(obj, None)
        if simplified_member:
            obj_name = simplified_member[0]
        else:
            return "No active bans found."
    title = f"{emojis.get_random_emoji()} Bans for {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary = await list_service.process_dictionary(
        cls=BanDictionary, dictionary=dictionary
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        ban_n = 0
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for member_snowflake, ban_dictionary in guild_data.get("members").items():
            member = guild.get_member(member_snowflake)
            if member:
                if not isinstance(obj, discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(url=obj.display_avatar.url)
                    thumbnail = True
            else:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake
                )
                if simplified_member:
                    display_name = simplified_member[0]
                    lines.append(f"**User:** {display_name} ({member_snowflake})")
                else:
                    continue
            for channel_snowflake, channel_dictionary in ban_dictionary.get(
                "bans"
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    continue
                if not isinstance(obj, discord.abc.GuildChannel):
                    lines.append(f"**Channel:** {channel.mention}")
                if isinstance(obj, discord.Member):
                    lines.append(
                        f"**Expires in:** {duration_builder.from_timestamp(channel_dictionary['expires_in']).to_unix_ts()}"
                    )
                    lines.append(f"**Reason:** {channel_dictionary['reason']}")
                    lines.append(
                        f"**Blacklisted:** {channel_dictionary['blacklisted']}"
                    )
                ban_n += 1
                field_count += 1
                if field_count >= list_service.CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = list_service.flush_page(embed, pages, title, guild.name)
                    lines = []
                    field_count = 0
        if lines:
            embed.add_field(name="Information", value="\n".join(lines), inline=False)
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({ban_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_guilds)
        pages.extend(processed_dictionary.skipped_members)
    if not pages:
        return "No bans found."
    return pages


async def build_blacklist_dictionary(obj) -> dict:
    database_factory = DatabaseFactory(MODEL)
    bans = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        bans = await database_factory.select(guild_snowflake=obj.id, singular=False)
    elif isinstance(obj, discord.abc.GuildChannel):
        bans = await database_factory.select(channel_snowflake=obj.id, singular=False)
    elif isinstance(obj, discord.Member):
        bans = await database_factory.select(member_snowflake=obj.id, singular=False)
    else:
        bans = await database_factory.select(singular=False)
    if bans:
        for ban in bans:
            if ban.blacklisted:
                dictionary.setdefault(ban.guild_snowflake, {"members": {}})
                dictionary[ban.guild_snowflake]["members"].setdefault(
                    ban.member_snowflake, {"bans": {}}
                )
                dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                    "bans"
                ].setdefault(ban.channel_snowflake, {})
    return dictionary


async def build_blacklist_pages(is_at_home: bool, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []
    thumbnail = False

    obj_name = "All Servers"
    if not isinstance(obj, int):
        obj_name = obj.name
    else:
        simplified_member = bot.registry.get(MemberState).active.get(obj, None)
        if simplified_member:
            obj_name = simplified_member[0]
        else:
            return "No active bans found."
    title = f"{emojis.get_random_emoji()} Blacklists for {obj_name}"

    dictionary = await build_blacklist_dictionary(obj=obj)
    processed_dictionary = await list_service.process_dictionary(
        cls=BanDictionary, dictionary=dictionary
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        ban_n = 0
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for member_snowflake, ban_dictionary in guild_data.get("members").items():
            member = guild.get_member(member_snowflake)
            if member:
                if not isinstance(obj, discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(url=obj.display_avatar.url)
                    thumbnail = True
            else:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if simplified_member:
                    display_name = simplified_member[0]
                    lines.append(f"**User:** {display_name} ({member_snowflake})")
                else:
                    continue
            for channel_snowflake, channel_dictionary in ban_dictionary.get(
                "bans"
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    continue
                if not isinstance(obj, discord.abc.GuildChannel):
                    lines.append(f"**Channel:** {channel.mention}")
                if isinstance(obj, discord.Member):
                    lines.append("**Blacklisted:** True")
                ban_n += 1
                field_count += 1
                if field_count >= list_service.CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = list_service.flush_page(embed, pages, title, guild.name)
                    lines = []
                    field_count = 0
        if lines:
            embed.add_field(name="Information", value="\n".join(lines), inline=False)
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({ban_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_guilds)
        pages.extend(processed_dictionary.skipped_members)
    if not pages:
        return "No blacklists found."
    return pages
