"""ban.py The purpose of this program is to inherit from DatabaseFactory to provide the ban moderation.

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

from datetime import datetime
from typing import Optional

import discord

from vyrtuous.fields.duration import DurationObject
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_guild_dictionary,
    flush_page,
)


class Ban(DatabaseFactory):

    ACT = "ban"
    CATEGORY = "ban"
    PLURAL = "Bans"
    SCOPES = ["channel", "member"]
    SINGULAR = "Ban"
    UNDO = "unban"

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "member_snowflake",
        "role_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "expires_in", "reason", "updated_at"]

    TABLE_NAME = "active_bans"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        role_snowflake: int,
        created_at: Optional[datetime] = None,
        expired: bool = False,
        expires_in: Optional[datetime] = None,
        reason: Optional[str] = "No reason provided.",
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.expired = expired
        self.expires_in = expires_in
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.reason = reason
        self.role_snowflake = role_snowflake
        self.updated_at = updated_at

    @classmethod
    async def act_embed(cls, action_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(action_information["action_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        role = source.guild.get_role(action_information["action_role_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been banned",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {action_information['action_duration']}\n"
                f"**Reason:** {action_information['action_reason']}\n"
                f"**Role:** {role.mention}"
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
        role = source.guild.get_role(action_information["action_role_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name}'s ban has been removed",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        kwargs = object_dict.get("columns", None)
        title = f"{get_random_emoji()} {Ban.PLURAL} {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        bans = await Ban.select(**kwargs)

        for ban in bans:
            guild_dictionary.setdefault(ban.guild_snowflake, {"members": {}})
            guild_dictionary[ban.guild_snowflake]["members"].setdefault(
                ban.member_snowflake, {"bans": {}}
            )
            guild_dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                "bans"
            ].setdefault(ban.channel_snowflake, {})
            guild_dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                "bans"
            ][ban.channel_snowflake] = {
                "reason": ban.reason,
                "expires_in": DurationObject.from_expires_in(ban.expires_in),
            }

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, ban_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in ban_dictionary.get(
                    "bans"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
                        lines.append(
                            f"**Expires in:** {channel_dictionary['expires_in']}"
                        )
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
                    field_count += 1
                    if field_count >= chunk_size:
                        embed.add_field(
                            name="Information", value="\n".join(lines), inline=False
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
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
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )
        return pages


class BanRole(DatabaseFactory):

    ACT = "brole"
    CATEGORY = "brole"
    PLURAL = "Ban Roles"
    SCOPES = ["channel", "guild"]
    SINGULAR = "Ban Role"
    UNDO = "brole"

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "role_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "updated_at"]

    TABLE_NAME = "ban_roles"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        role_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.created_at = created_at
        self.channel_snowflake = channel_snowflake
        self.guild_snowflake = guild_snowflake
        self.role_snowflake = role_snowflake
        self.updated_at = updated_at
