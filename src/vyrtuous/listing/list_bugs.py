"""!/bin/python3
bug_service.py The purpose of this program is to extend Service to service the bug command class.

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
from vyrtuous.db.bug import Bug
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = Bug


@dataclass
class BugDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Any]]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_messages: List[discord.Embed] = field(default_factory=list)


async def build_dictionary(obj, reference):
    database_factory = DatabaseFactory(MODEL)
    bugs = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        bugs = await database_factory.select(
            guild_snowflake=obj.id, reference=reference, singular=False
        )
    elif not obj and reference:
        bugs = await database_factory.select(reference=reference, singular=False)
    else:
        bugs = await database_factory.select(singular=False)
    if bugs:
        for bug in bugs:
            dictionary.setdefault(bug.guild_snowflake, {"messages": {}})
            messages = dictionary[bug.guild_snowflake]["messages"]
            messages.setdefault(
                bug.message_snowflake,
                {
                    "channel_snowflake": bug.channel_snowflake,
                    "developer_snowflakes": [],
                    "id": bug.id,
                    "notes": [],
                    "resolved": bug.resolved,
                },
            )
            messages[bug.message_snowflake]["developer_snowflakes"].extend(
                bug.member_snowflakes
            )
            messages[bug.message_snowflake]["notes"].append(bug.notes)
    return dictionary


async def build_pages(is_at_home: bool, obj, reference, scope):
    bot = DiscordBot.get_instance()
    lines, pages = [], []

    obj_name = "All Servers"
    if obj:
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Developer Logs in {obj_name}"

    dictionary = await build_dictionary(obj=obj, reference=reference)
    processed_dictionary = await list_service.process_dictionary(
        cls=BugDictionary, dictionary=dictionary
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        bug_n = 0
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for message_snowflake, entry in guild_data.get("messages", {}).items():
            channel = guild.get_channel(entry["channel_snowflake"])
            if channel is None:
                continue
            if not isinstance(channel, discord.TextChannel) or not isinstance(
                channel, discord.VoiceChannel
            ):
                continue
            if scope == "resolved" and not entry.get("resolved"):
                continue
            if scope == "unresolved" and entry.get("resolved"):
                continue
            msg = await channel.fetch_message(message_snowflake)
            lines.append(
                f"**Resolved:** {'\u2705' if entry.get('resolved') else '\u274c'}"
            )
            lines.append(f"**Message:** {msg.jump_url}")
            if reference == str(entry["id"]):
                lines.append(
                    f"**Notes:** {entry['notes'] if entry.get('notes') is not None else None}"
                )
                lines.append(
                    f"**Assigned to:** {', '.join(str(d) for d in entry['developer_snowflakes']) if entry.get('developer_snowflakes') else None}"
                )
            else:
                lines.append(f"**Reference:** {entry['id']}")
            bug_n += 1
            field_count += 1
            if field_count >= list_service.CHUNK_SIZE:
                embed.add_field(
                    name=f"**Channel:** {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed = list_service.flush_page(embed, pages, title, guild.name)
                lines = []
            if lines:
                embed.add_field(
                    name=f"**Channel:** {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({bug_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_guilds)
        pages.extend(processed_dictionary.skipped_messages)
    if not pages:
        return "No bugs found."
    return pages
