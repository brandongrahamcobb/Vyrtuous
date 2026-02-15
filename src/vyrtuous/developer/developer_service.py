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

from vyrtuous.developer.developer import Developer
from vyrtuous.sysadmin.sysadmin_service import SysadminService


class NotDeveloper(commands.CommandError):
    def __init__(self, message="Member is not a developer."):
        super().__init__(message)


class DeveloperService:
    __CHUNK_SIZE = 7
    MODEL = Developer

    def __init__(
        self, *, author_service=None, bot=None, database_factory=None, emoji=None
    ):
        self.__author_service = author_service
        self.__bot = bot
        self.__database_factory = database_factory
        self.__database_factory = self.MODEL
        self.__emoji = emoji
        self.__sysadmin_service = SysadminService(
            author_service=author_service, bot=bot
        )

    async def is_developer(self, member_snowflake: int) -> bool:
        developer = await self.__database_factory.select(
            member_snowflake=int(member_snowflake), singular=True
        )
        if not developer:
            raise NotDeveloper
        return True

    async def is_developer_wrapper(
        self, source: Union[commands.Context, discord.Interaction, discord.Message]
    ):
        member = self.__author_service.resolve_author(source=source)
        member_snowflake = member.id
        return await self.is_developer(member_snowflake)

    async def build_clean_dictionary(self, where_kwargs):
        dictionary = {}
        developers = await self.__database_factory.select(
            singular=False, **where_kwargs
        )
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

    async def build_pages(self, object_dict, **kwargs):
        lines = []
        pages = []
        title = f"{self.__emoji.get_random_emoji()} Developers {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), (discord.Guild, discord.Member)) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(where_kwargs=where_kwargs)

        embed = discord.Embed(
            title=title, description="All guilds", color=discord.Color.blue()
        )
        dev_n = 0
        for key, values in dictionary.items():
            field_count = 0
            thumbnail = False
            for member_snowflake, member_data in values.items():
                user = self.__bot.get_user(member_snowflake)
                if not user:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {user.display_name} {user.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                dev_n += 1
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    pages.append(embed)
                    embed = (
                        discord.Embed(
                            title=title,
                            description="All guilds continued...",
                            color=discord.Color.blue(),
                        ),
                    )
                    field_count = 0
                    lines = []
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)
        if pages:
            pages[0].description = f"All guilds **({dev_n})**"
        return pages

    async def report_issue(self, message, reference, source, user):
        online_developer_mentions = []
        member = self.__bot.get_user(self.__config["discord_owner_id"])
        online_developer_mentions.append(member.mention)
        if source.guild:
            developers = await self.__database_factory.select(
                guild_snowflake=source.guild.id
            )
            message = f"Issue reported by {user.name}!\n**Message:** {message.jump_url}\n**Reference:** {reference}"
            for dev in developers:
                member = source.guild.get_member(dev.member_snowflake)
                if member and member.status != discord.Status.offline:
                    online_developer_mentions.append(member.mention)
                    try:
                        await member.send(message)
                    except discord.Forbidden as e:
                        self.__bot.logger.warning(
                            f"Unable to send a developer log ID: {id}. {str(e).capitalize()}"
                        )
        message = "Your report has been submitted"
        if online_developer_mentions:
            message = f"{message}. The developers {', '.join(online_developer_mentions)} are online and will respond to your report shortly."
        await user.send(message)

    async def toggle_developer(self, member_dict):
        where_kwargs = member_dict.get("columns", None)
        del where_kwargs["guild_snowflake"]
        developer = await self.__database_factory.select(singular=True, **where_kwargs)
        if developer:
            await self.__database_factory.delete(**where_kwargs)
            action = "revoked"
        else:
            developer = self.MODEL(**where_kwargs)
            await self.__database_factory.create(developer)
            action = "granted"
        return f"Developer access for {member_dict.get('mention', None)} has been {action} globally."
