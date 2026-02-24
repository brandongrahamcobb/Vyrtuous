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

from copy import copy

import discord
from discord.ext import commands

from vyrtuous.developer.developer import Developer


class NotDeveloper(commands.CheckFailure):
    def __init__(self, message="Member is not a developer."):
        super().__init__(message)


class DeveloperService:
    __CHUNK_SIZE = 12
    MODEL = Developer

    def __init__(
        self,
        *,
        bot=None,
        bug_service=None,
        database_factory=None,
        duration_builder=None,
        emoji=None,
    ):
        self.__bot = bot
        self.__bug_service = bug_service
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__duration_builder = duration_builder
        self.__emoji = emoji

    async def is_developer(self, member_snowflake: int) -> bool:
        developer = await self.__database_factory.select(
            member_snowflake=int(member_snowflake), singular=True
        )
        if not developer:
            raise NotDeveloper
        return True

    async def is_developer_wrapper(self, context):
        return await self.is_developer(member_snowflake=int(context.author.id))

    async def build_dictionary(self, kwargs):
        dictionary = {}
        developers = await self.__database_factory.select(singular=False, **kwargs)
        for developer in developers:
            dictionary.setdefault("members", {})
            dictionary["members"].setdefault(
                developer.member_snowflake, {"developers": {}}
            )
            dictionary["members"][developer.member_snowflake]["developers"].update(
                {"placeholder": "placeholder"}
            )
        return dictionary

    async def build_pages(self, object_dict):
        lines, pages = [], []

        obj = object_dict.get("object")
        obj_name = "All Servers"
        if isinstance(obj, discord.Guild):
            obj_name = obj.name
        elif isinstance(obj, discord.Member):
            obj_name = object_dict.get("name", None)
        title = f"{self.__emoji.get_random_emoji()} Developers for {obj_name}"

        dictionary = await self.build_dictionary(
            kwargs=object_dict.get("columns", None)
        )

        embed = discord.Embed(
            title=title, description="Information", color=discord.Color.blue()
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
                            description=title,
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
            original_description = embed.description or ""
            embed.description = f"**{original_description}** **({dev_n})**"
            pages.append(embed)
        return pages

    async def report_issue(self, author, message, reference, source):
        online_developer_mentions = []
        sysadmin = self.__bot.get_user(self.__bot.config["discord_owner_id"])
        online_developer_mentions.append(sysadmin.mention)
        if source.guild:
            message = f"Issue reported by {author.name}!\n**Message:** {message.jump_url}\n**Reference:** {reference}"
            for developer in await self.developers():
                member = source.guild.get_member(developer.member_snowflake)
                if member and member.status != discord.Status.offline:
                    online_developer_mentions.append(member.mention)
                    try:
                        await member.send(message)
                        await sysadmin.send(message)
                    except discord.Forbidden as e:
                        self.__bot.logger.warning(
                            f"Unable to send a developer log ID: {id}. {str(e).capitalize()}"
                        )
        message = "Your report has been submitted"
        if online_developer_mentions:
            message = f"{message}. The developers {', '.join(online_developer_mentions)} are online and will respond to your report shortly."
        await author.send(message)

    async def toggle_developer(self, member_dict):
        kwargs = member_dict.get("columns", None)
        del kwargs["guild_snowflake"]
        found = False
        for developer in await self.developers():
            if developer.member_snowflake == member_dict.get("id", None):
                await self.__database_factory.delete(**kwargs)
                found = True
                action = "revoked"
        if not found:
            new_developer = self.MODEL(**kwargs)
            await self.__database_factory.create(new_developer)
            action = "granted"
        return f"Developer access for {member_dict.get('mention', None)} has been {action} globally."

    async def ping_about_expired_bugs(
        self,
        channel,
        embed,
        member,
        member_snowflakes,
        msg,
        notes,
        updated_at,
    ):
        assigned_developer_mentions = []
        for developer_snowflake in member_snowflakes:
            assigned_developer = self.__bot.get_user(developer_snowflake)
            assigned_developer_mentions.append(assigned_developer.mention)
        embed.add_field(
            name=f"Updated: {self.__duration_builder.from_timestamp(updated_at).to_unix_ts()}",
            value=f"**Link:** {msg.jump_url}\n**Developers:** {', '.join(assigned_developer_mentions)}\n**Notes:** {notes}",
            inline=False,
        )
        for developer_snowflake in member_snowflakes:
            member = self.__bot.get_user(developer_snowflake)
            if member is None:
                self.__bot.logger.info(
                    f"Unable to locate member {member.id} in channel {channel.name} ({channel.id}) not sending developer log."
                )
                continue
            try:
                await member.send(embed=embed)
                self.__bot.logger.info(
                    f"Sent the issue to member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id})."
                )
            except Exception as e:
                self.__bot.logger.warning(
                    f"Unable to send the issue to member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}). {str(e).capitalize()}"
                )

    async def developers(self):
        developers = await self.__database_factory.select()
        return developers

    async def handle_developer_assignment(self, member_dict, reference):
        for developer in await self.developers():
            if developer.member_snowflake == member_dict.get("id", None):
                bug, state = self.__bug_service.handle_bug_assignment(
                    developer=developer, reference=reference
                )
                if not state:
                    action = "unassigned"
                else:
                    action = "assigned"
                return await self.__bug_service.create_embed(
                    action=action,
                    bug=bug,
                    member_dict=member_dict,
                )
