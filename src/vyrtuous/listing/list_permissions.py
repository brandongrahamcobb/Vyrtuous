"""!/bin/python3
permission_service.py The purpose of this program is to provide the service for deciding whether a member has sufficient permissions.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

import discord

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


def build_dictionary(obj) -> dict[int, dict[str, dict[int, dict[str, list[str]]]]]:
    channels = []
    dictionary: dict[int, dict[str, dict[int, dict[str, list[str]]]]] = {}
    if isinstance(
        obj, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
    ):
        me = obj.guild.me
        channels = [obj]
    else:
        me = obj.me
        channels = [channel for channel in obj.channels]
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


async def build_pages(
    obj,
) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = obj.name
    if bot.user is None:
        return "The bot must be in a server."
    title = f"{emojis.get_random_emoji()} {bot.user.display_name} Missing Permissions in {obj_name}"

    dictionary = build_dictionary(obj=obj)
    processed_dictionary: list_service.PermissionDictionary = (
        await list_service.process_dictionary(
            cls=list_service.PermissionDictionary, dictionary=dictionary
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
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
    return pages
