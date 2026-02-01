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


import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.roles.developer import Developer
from vyrtuous.utils.logger import logger
from vyrtuous.utils.dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_messages,
    clean_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.inc.helpers import CHUNK_SIZE


class Bug(DatabaseFactory):

    ACT = "bug"
    CATEGORY = "bug"
    PLURAL = "Bugs"
    SCOPES = ["channels"]
    SINGULAR = "Bug"
    UNDO = None

    REQUIRED_ARGS = [
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
    lines, pages = [], []

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        id: str,
        member_snowflakes: list[int],
        message_snowflake: int,
        notes: str | None,
        resolved: bool = False,
        created_at: datetime = datetime.now(timezone.utc),
        updated_at: datetime = datetime.now(timezone.utc),
    ):
        self.channel_snowflake = channel_snowflake
        self.created_at: datetime = created_at or datetime.now(timezone.utc)
        self.member_snowflakes = member_snowflakes
        self.id = id
        self.guild_snowflake = guild_snowflake
        self.message_snowflake = message_snowflake
        self.notes = notes
        self.resolved = resolved
        self.updated_at: datetime = updated_at or datetime.now(timezone.utc)

    async def create_embed(self, action, member_snowflake):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(self.channel_snowflake)
        guild = bot.get_guild(self.guild_snowflake)
        try:
            msg = await channel.fetch_message(self.message_snowflake)
        except discord.NotFound:
            logger.warning(f"Message reference not found ({self.message_snowflake}).")
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
                f"**Assigned devs:** {', '.join(user_mentions)}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        bugs = await Bug.select(singular=False, **where_kwargs)
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
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_messages = await generate_skipped_messages(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_messages=skipped_messages,
        )
        if is_at_home:
            if skipped_guilds:
                Bug.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_messages:
                Bug.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_messages,
                        title="Skipped Messages in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, filter, where_kwargs, is_at_home):
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Developer Logs"

        dictionary = await Bug.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
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
                msg = await channel.fetch_message(message_snowflake)
                Bug.lines.append(
                    f'**Resolved:** {"\u2705" if entry.get("resolved") else "\u274c"}'
                )
                Bug.lines.append(f"**Message:** {msg.jump_url}")
                if where_kwargs.get("id", None) == str(entry["id"]):
                    Bug.lines.append(
                        f'**Notes:** {entry["notes"] if entry.get("notes") is not None else None}'
                    )
                    Bug.lines.append(
                        f'**Assigned to:** {", ".join(str(d) for d in entry["developer_snowflakes"]) if entry.get("developer_snowflakes") else None}'
                    )
                else:
                    Bug.lines.append(f'**Reference:** {entry["id"]}')
                field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name=f"**Channel:** {channel.mention}",
                        value="\n".join(Bug.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, Bug.pages, title, guild.name)
                    Bug.lines = []
            if Bug.lines:
                embed.add_field(
                    name=f"**Channel:** {channel.mention}",
                    value="\n".join(Bug.lines),
                    inline=False,
                )
            Bug.pages.append(embed)
        return Bug.pages

    @classmethod
    async def assign_bug_to_developer(cls, reference, member_dict):
        where_kwargs = member_dict.get("columns", None)

        developer = await Developer.select(singular=True, **where_kwargs)
        if not developer:
            return (
                f"Developer not found for target ({member_dict.get("mention", None)})."
            )

        bug = await Bug.select(id=reference, resolved=False, singular=True)
        if not bug:
            return f"Unresolved issue not found for reference: {reference}."
        member_snowflakes = bug.member_snowflakes
        where_kwargs = {"id": bug.id}
        member_snowflakes = bug.member_snowflakes
        if developer.member_snowflake in bug.member_snowflakes:
            member_snowflakes.remove(developer.member_snowflake)
            set_kwargs = {"member_snowflakes": member_snowflakes}
            await bug.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            embed = await bug.create_embed(
                action="unassigned",
                member_snowflake=developer.member_snowflake,
            )
            return embed
        else:
            member_snowflakes.append(developer.member_snowflake)
            set_kwargs = {"member_snowflakes": member_snowflakes}
            await bug.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            embed = await bug.create_embed(
                action="assigned",
                member_snowflake=developer.member_snowflake,
            )
            await member_dict.get("object", None).send(embed=embed)
            return embed
