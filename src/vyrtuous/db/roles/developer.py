"""developer.py The purpose of this program is to inherit from the DatabaseFactory to provide the developer role.

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
from typing import Union

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.roles.sysadmin import is_sysadmin_wrapper
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.dir_to_classes import skip_db_discovery
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.inc.helpers import CHUNK_SIZE


@skip_db_discovery
class NotDeveloper(commands.CheckFailure):
    def __init__(self, message="Member is not a sysadmin or developer."):
        super().__init__(message)


async def is_developer_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_developer(member_snowflake)


def developer_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (is_sysadmin_wrapper, is_developer_wrapper):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure("Member is not a sysadmin or developer.")

    predicate._permission_level = "Developer"
    return commands.check(predicate)


async def is_developer(member_snowflake: int) -> bool:
    developer = await Developer.select(member_snowflake=int(member_snowflake), singular=True)
    if not developer:
        raise NotDeveloper
    return True


class Developer(DatabaseFactory):

    ACT = "dev"
    CATEGORY = "dev"
    PLURAL = "Developers"
    SCOPES = ["member"]
    SINGULAR = "Developer"
    UNDO = "dev"

    REQUIRED_ARGS = ["member_snowflake"]
    OPTIONAL_ARGS = ["created_at", "updated_at"]

    TABLE_NAME = "developers"
    lines, pages = [], []

    def __init__(
        self,
        member_snowflake: int,
        created_at: datetime = datetime.now(timezone.utc),
        updated_at: datetime = datetime.now(timezone.utc),
        **kwargs,
    ):
        self.created_at = created_at
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.updated_at = updated_at

    @classmethod
    async def build_clean_dictionary(cls, where_kwargs):
        dictionary = {}
        developers = await Developer.select(singular=False, **where_kwargs)
        for developer in developers:
            dictionary.setdefault("members", {})
            dictionary["members"].setdefault(
                developer.member_snowflake, {"developers": {}}
            )
            dictionary["members"][developer.member_snowflake]["developers"].update(
                {"placeholder": "placeholder"}
            )
        cleaned_dictionary = dictionary
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, **kwargs):
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} {Developer.PLURAL} {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), (discord.Guild, discord.Member)) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await Developer.build_clean_dictionary(where_kwargs=where_kwargs)

        embed = discord.Embed(
            title=title, description="All guilds", color=discord.Color.blue()
        )
        for key, values in dictionary.items():
            field_count = 0
            thumbnail = False
            for member_snowflake, member_data in values.items():
                user = bot.get_user(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    Developer.lines.append(
                        f"**User:** {user.display_name} {user.mention}"
                    )
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(Developer.lines),
                        inline=False,
                    )
                    Developer.pages.append(embed)
                    embed = (
                        discord.Embed(
                            title=title,
                            description="All guilds continued...",
                            color=discord.Color.blue(),
                        ),
                    )
                    field_count = 0
                    Developer.lines = []
            if Developer.lines:
                embed.add_field(
                    name="Information", value="\n".join(Developer.lines), inline=False
                )
            Developer.pages.append(embed)
        return Developer.pages

    @classmethod
    async def toggle_developer(cls, member_dict, snowflake_kwargs):
        bot = DiscordBot.get_instance()
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        guild = bot.get_guild(guild_snowflake)
        where_kwargs = member_dict.get("columns", None)
        developer = await Developer.select(singular=True, **where_kwargs)
        if developer:
            await Developer.delete(**where_kwargs)
            action = "revoked"
        else:
            developer = Developer(**where_kwargs)
            await developer.create()
            action = "granted"
        return f"Developer access for {member_dict.get("mention", None)} has been {action} in {guild.name}."
