"""flag.py The purpose of this program is to inherit from DatabaseFactory to provide the flag moderation.

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
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_dictionary,
    flush_page,
)
from vyrtuous.inc.helpers import CHUNK_SIZE


class Flag(DatabaseFactory):

    ACT = "flag"
    CATEGORY = "flag"
    PLURAL = "Flags"
    SCOPES = ["channel", "member"]
    SINGULAR = "Flag"
    UNDO = "unflag"

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "member_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "reason", "updated_at"]

    TABLE_NAME = "active_flags"
    lines, pages = [], []

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: datetime = datetime.now(timezone.utc),
        reason: str = "No reason provided.",
        updated_at: datetime = datetime.now(timezone.utc),
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.reason = reason
        self.updated_at = updated_at

    @classmethod
    async def act_embed(cls, action_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(action_information["action_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been flagged",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Reason:** {action_information['action_reason']}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, action_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(action_information["action_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name}'s flag has been removed",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        flags = await Flag.select(**where_kwargs)
        for flag in flags:
            dictionary.setdefault(flag.guild_snowflake, {"members": {}})
            dictionary[flag.guild_snowflake]["members"].setdefault(
                flag.member_snowflake, {"flags": {}}
            )
            dictionary[flag.guild_snowflake]["members"][flag.member_snowflake][
                "flags"
            ].setdefault(flag.channel_snowflake, {})
            dictionary[flag.guild_snowflake]["members"][flag.member_snowflake]["flags"][
                flag.channel_snowflake
            ].update(
                {
                    "reason": flag.reason,
                }
            )
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                Flag.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                Flag.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} {Flag.PLURAL} {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await Flag.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, flag_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    Flag.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in flag_dictionary.get(
                    "flags", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        Flag.lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
                        Flag.lines.append(f"**Reason:** {channel_dictionary['reason']}")
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(Flag.lines),
                            inline=False,
                        )
                        embed = flush_page(embed, Flag.pages, title, guild.name)
                        Flag.lines = []
            if Flag.lines:
                embed.add_field(
                    name="Information", value="\n".join(Flag.lines), inline=False
                )
            Flag.pages.append(embed)
        return Flag.pages
