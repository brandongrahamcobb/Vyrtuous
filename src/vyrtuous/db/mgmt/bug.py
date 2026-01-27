"""developer.py The purpose of this program is to inherit from the user class as a developer.

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
from uuid import UUID
from typing import Optional

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.logger import logger
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_messages,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji

class Bug(DatabaseFactory):

    ACT = "bug"
    CATEGORY = "bug"
    PLURAL = "Bugs"
    SCOPES = ["channels"]
    SINGULAR = "Bug"
    UNDO = None

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "id",
        "member_snowflakes",
        "message_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "notes",
        "resolved",
        "updated_at",
    ]
    TABLE_NAME = "bug_tracking"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        id: str,
        member_snowflakes: list[int],
        message_snowflake: int,
        created_at: Optional[datetime] = None,
        notes: Optional[str] = None,
        resolved: Optional[bool] = False,
        updated_at: Optional[datetime] = None,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.created_at: datetime = created_at or datetime.now(timezone.utc)
        self.member_snowflakes = member_snowflakes
        self.id = id
        self.guild_snowflake = guild_snowflake
        self.message_snowflake = message_snowflake
        self.notes = notes
        self.resolved = resolved
        self.updated_at: datetime = updated_at or datetime.now(timezone.utc)

    async def create_embed(self, action, member_snowflake, source):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(self.channel_snowflake)
        guild = bot.get_guild(self.guild_snowflake)
        try:
            msg = await channel.fetch_message(self.message_snowflake)
        except discord.NotFound:
            logger.warning(f"Message reference not found ({self.message_snowflake}).")
        author = resolve_author(source=source)
        member = bot.get_user(member_snowflake)
        user_mentions = []
        for member_snowflake in self.member_snowflakes:
            user_mentions.append(bot.get_user(member_snowflake).mention)
        embed = discord.Embed(
            title=f"{get_random_emoji()} {member.display_name} has been {action}",
            description=(
                f"**Guild:** {guild.name}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Message:** {msg.jump_url}\n"
                f"**By:** {author.mention}\n"
                f"**Assigned devs:** {', '.join(user_mentions)}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
    
    @classmethod
    async def build_pages(cls, filter, kwargs, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} Developer Logs"

        bugs = await Bug.select(**kwargs)

        for bug in bugs:
            guild_dictionary.setdefault(bug.guild_snowflake, {"messages": {}})
            messages = guild_dictionary[bug.guild_snowflake]["messages"]
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

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_messages = await generate_skipped_messages(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_messages=skipped_messages,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for message_snowflake, entry in guild_data.get("messages", {}).items():
                channel = bot.get_channel(entry["channel_snowflake"])
                if not channel:
                    continue
                if filter == "resolved" and not entry.get("resolved"):
                    continue
                if filter == "unresolved" and entry.get("resolved"):
                    continue
                lines = []
                msg = await channel.fetch_message(message_snowflake)
                lines.append(
                    f'**Resolved:** {"\u2705" if entry.get("resolved") else "\u274c"}'
                )
                lines.append(f"**Message:** {msg.jump_url}")
                if kwargs.get("id", None) == str(entry["id"]):
                    lines.append(
                        f'**Notes:** {entry["notes"] if entry.get("notes") is not None else None}'
                    )
                    lines.append(
                        f'**Assigned to:** {", ".join(str(d) for d in entry["developer_snowflakes"]) if entry.get("developer_snowflakes") else None}'
                    )
                else:
                    lines.append(f'**Reference:** {entry["id"]}')
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name=f"**Channel:** {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name=f"**Channel:** {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_messages:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_messages,
                    title="Skipped Messages in Server",
                )
        return pages