"""!/bin/python3
developer_service.py The purpose of this program is to extend Service to service the developer class.

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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.base.record_service import RecordService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.developer.developer import Developer
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.sysadmin.sysadmin_service import is_sysadmin_wrapper
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.errors import NotDeveloper


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
        raise NotDeveloper

    predicate._permission_level = "Developer"
    return commands.check(predicate)


async def is_developer(member_snowflake: int) -> bool:
    developer = await Developer.select(
        member_snowflake=int(member_snowflake), singular=True
    )
    if not developer:
        raise NotDeveloper
    return True


class DeveloperService(RecordService):
    lines, pages = [], []
    model = Developer

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
        title = f"{get_random_emoji()} Developers {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), (discord.Guild, discord.Member)) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await DeveloperService.build_clean_dictionary(
            where_kwargs=where_kwargs
        )

        embed = discord.Embed(
            title=title, description="All guilds", color=discord.Color.blue()
        )
        dev_n = 0
        for key, values in dictionary.items():
            field_count = 0
            thumbnail = False
            for member_snowflake, member_data in values.items():
                user = bot.get_user(member_snowflake)
                if not user:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    DeveloperService.lines.append(
                        f"**User:** {user.display_name} {user.mention}"
                    )
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                dev_n += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(DeveloperService.lines),
                        inline=False,
                    )
                    DeveloperService.pages.append(embed)
                    embed = (
                        discord.Embed(
                            title=title,
                            description="All guilds continued...",
                            color=discord.Color.blue(),
                        ),
                    )
                    field_count = 0
                    DeveloperService.lines = []
            if DeveloperService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(DeveloperService.lines),
                    inline=False,
                )
            DeveloperService.pages.append(embed)
        if DeveloperService.pages:
            DeveloperService.pages[0].description = f"All guilds **({dev_n})**"
        return DeveloperService.pages

    @classmethod
    async def toggle_developer(cls, member_dict):
        where_kwargs = member_dict.get("columns", None)
        del where_kwargs["guild_snowflake"]
        developer = await Developer.select(singular=True, **where_kwargs)
        if developer:
            await Developer.delete(**where_kwargs)
            action = "revoked"
        else:
            developer = Developer(**where_kwargs)
            await developer.create()
            action = "granted"
        return f"Developer access for {member_dict.get('mention', None)} has been {action} globally."
