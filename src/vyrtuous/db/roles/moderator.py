"""moderator.py The purpose of this program is to inherit from the DatabaseFactory to provide the moderator role.

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
from typing import Optional, Union

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.roles.administrator import is_administrator_wrapper
from vyrtuous.db.roles.coordinator import is_coordinator_wrapper
from vyrtuous.db.roles.developer import is_developer_wrapper
from vyrtuous.db.roles.guild_owner import is_guild_owner_wrapper
from vyrtuous.db.roles.sysadmin import is_sysadmin_wrapper
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.dir_to_classes import skip_db_discovery
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji


@skip_db_discovery
class NotModerator(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a sysadmin, developer, guild owner, administrator, coordinator or moderator in this channel.",
    ):
        super().__init__(message)


async def is_moderator_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_moderator(
        channel_snowflake=source.channel.id,
        guild_snowflake=source.guild.id,
        member_snowflake=member_snowflake,
    )


async def is_moderator(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    moderator = await Moderator.select(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        singular=True,
    )
    if not moderator:
        raise NotModerator
    return True


def moderator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_sysadmin_wrapper,
            is_developer_wrapper,
            is_guild_owner_wrapper,
            is_administrator_wrapper,
            is_coordinator_wrapper,
            is_moderator_wrapper,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "Member is not a sysadmin, developer, guild owner, administrator, coordinator or moderator in this channel."
        )

    predicate._permission_level = "Moderator"
    return commands.check(predicate)


class Moderator(DatabaseFactory):

    ACT = "mod"
    CATEGORY = "mod"
    PLURAL = "Moderators"
    SCOPES = ["channel", "member"]
    SINGULAR = "Moderator"
    UNDO = "mod"

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "member_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "updated_at"]

    TABLE_NAME = "moderators"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.updated_at = updated_at

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        title = f"{get_random_emoji()} {Moderator.PLURAL} {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"
        kwargs = object_dict.get("columns", None)

        moderators = await Moderator.select(**kwargs)

        for moderator in moderators:
            guild_dictionary.setdefault(moderator.guild_snowflake, {"members": {}})
            guild_dictionary[moderator.guild_snowflake]["members"].setdefault(
                moderator.member_snowflake, {"moderators": {}}
            )
            guild_dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"].setdefault(moderator.channel_snowflake, {})
            guild_dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"][moderator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )

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
            for member_snowflake, member_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in member_dictionary.get(
                    "moderators", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    lines.append(f"**Channel:** {channel.mention}")
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

    @classmethod
    async def toggle_moderator(cls, channel_dict, member_dict, snowflake_kwargs):
        from vyrtuous.utils.check import has_equal_or_lower_role

        await has_equal_or_lower_role(
            snowflake_kwargs=snowflake_kwargs,
            member_snowflake=member_dict.get("id", None),
        )
        where_kwargs = {}
        where_kwargs.update(channel_dict.get("columns", None))
        where_kwargs.update(member_dict.get("columns", None))
        moderator = await Moderator.select(**where_kwargs, singular=True)
        if moderator:
            await Moderator.delete(**where_kwargs)
            action = "revoked"
        else:
            moderator = Moderator(**where_kwargs)
            await moderator.create()
            action = "granted"
        return (
            f"Moderator access for {member_dict.get("mention", None)} has been "
            f"{action} in {channel_dict.get("mention", None)}."
        )
