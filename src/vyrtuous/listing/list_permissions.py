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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

TARGET_PERMISSIONS = (
    "add_reactions",
    "manage_messages",
    "move_members",
    "mute_members",
    "send_messages",
    "view_channel",
)


@dataclass(frozen=True)
class PermissionDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, List[str]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


def build_dictionary(obj, me) -> dict:
    bot: DiscordBot = DiscordBot.get_instance()
    channels = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        channels = obj.channels
    elif isinstance(obj, discord.abc.GuildChannel):
        channels = [obj]
    else:
        channels = [channel for guild in bot.guilds for channel in guild.channels]
    if channels:
        for channel in channels:
            permissions = channel.permissions_for(me)
            missing = []
            for permission in TARGET_PERMISSIONS:
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


async def build_pages(obj, context, is_at_home) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    guild = bot.get_guild(context.guild.id)
    if guild is None:
        return "This command must be used in a server."
    obj_name = "All Servers"
    if obj:
        obj_name = obj.name
    if bot.user is None:
        return "The bot must be in a server."
    title = f"{emojis.get_random_emoji()} {bot.user.display_name} Missing Permissions in {obj_name}"

    dictionary = build_dictionary(obj=obj, me=guild.me)
    processed_dictionary = await list_service.process_dictionary(
        cls=PermissionDictionary, dictionary=dictionary
    )

    for guild_snowflake, guild_data in dictionary.items():
        perm_n = 0
        field_count = 0
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for channel_snowflake, channel_data in guild_data.get("channels", {}).items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            lines.append(f"Channel: {channel.mention}")
            for _, permissions in channel_data.items():
                for permission in permissions:
                    lines.append(f"  ↳ {permission}")
            field_count += 1
            perm_n += 1
            if field_count >= list_service.CHUNK_SIZE:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
                embed = list_service.flush_page(embed, pages, title, guild.name)
                lines = []
        if lines:
            embed.add_field(
                name="Information",
                value="\n".join(lines),
                inline=False,
            )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({perm_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_channels)
        pages.extend(processed_dictionary.skipped_guilds)
    return pages
