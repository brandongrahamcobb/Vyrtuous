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

from copy import copy
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List

import discord

from vyrtuous.bug.bug import Bug


@dataclass
class BugDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Any]]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_messages: List[discord.Embed] = field(default_factory=list)


class BugService:
    __CHUNK_SIZE = 7
    MODEL = Bug

    def __init__(
        self,
        bot=None,
        database_factory=None,
        developer_service=None,
        dictionary_service=None,
        emoji=None,
    ):
        self.__bot = bot
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__developer_service = developer_service

    async def interact_with_bug(self, action, notes, reference):
        message = "You successfully "
        bug = await self.__database_factory.select(
            id=reference, resolved=False, singular=True
        )
        if not bug:
            return f"Unresolved issue not found for reference ({reference})."
        if action and action.lower() == "resolve":
            where_kwargs = {"id": reference}
            set_kwargs = {"resolved": True}
            await self.__database_factory.update(
                where_kwargs=where_kwargs, set_kwargs=set_kwargs
            )
            message += "resolved the issue. The record will remain in the database for the next 30 days."
        elif action and action.lower() == "append":
            where_kwargs = {"id": reference}
            set_kwargs = {"notes": bug.notes + notes if bug.notes else notes}
            await self.__database_factory.update(
                where_kwargs=where_kwargs, set_kwargs=set_kwargs
            )
            message += "appended to the previous notes."
        elif action and action.lower() == "overwrite":
            where_kwargs = {"id": reference}
            set_kwargs = {"notes": notes}
            await self.__database_factory.update(
                where_kwargs=where_kwargs, set_kwargs=set_kwargs
            )
            message += "overwrote the previous notes."
        return message

    async def clean_expired(self):
        now = datetime.now(timezone.utc)
        bugs = await self.__database_factory.select(resolved=True)
        if bugs:
            for bug in bugs:
                channel_snowflake = int(bug.channel_snowflake)
                guild_snowflake = int(bug.guild_snowflake)
                member_snowflakes = int(bug.member_snowflakes)
                message_snowflake = int(bug.message_snowflake)
                reference = bug.id
                if bug.created_at < now - timedelta(weeks=1):
                    guild = self.__bot.get_guild(bug.guild_snowflake)
                    if guild is None:
                        self.__bot.logger.info(
                            f"Unable to locate guild {guild_snowflake}, not sending developer log."
                        )
                        continue
                    embed = discord.Embed(
                        title=f"\U000026a0\U0000fe0f An issue is unresolved in {guild.name}",
                        color=discord.Color.red(),
                    )
                    channel = guild.get_channel(channel_snowflake)
                    if channel is None:
                        self.__bot.logger.info(
                            f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}, not sending developer log."
                        )
                        continue
                    for member_snowflake in member_snowflakes:
                        member = channel.get_member(member_snowflake)
                        if member is None:
                            self.__bot.logger.info(
                                f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), not sending developer log."
                            )
                            continue
                        embed.set_thumbnail(url=member.display_avatar.url)
                        try:
                            msg = await channel.fetch_message(message_snowflake)
                        except Exception as e:
                            self.__bot.logger.warning(
                                f"Unable to locate a message {msg} in {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), deleting developer log. {str(e).capitalize()}"
                            )
                            return await self.__database_factory.delete(id=reference)
                        time_since_updated = await self.duration.from_expires_at(
                            bug.updated_at
                        )
                        self.__developer_service.clean_expired(
                            channel=channel,
                            embed=embed,
                            guild_snowflake=guild_snowflake,
                            member=member,
                            member_snowflakes=bug.member_snowflakes,
                            msg=msg,
                            notes=bug.notes,
                            time_since_updated=time_since_updated,
                        )

    async def create_embed(self, action, bug, member_snowflake):
        guild = self.__bot.get_guild(bug.guild_snowflake)
        channel = guild.get_channel(bug.channel_snowflake)
        try:
            msg = await channel.fetch_message(bug.message_snowflake)
        except discord.NotFound:
            self.__bot.logger.warning(
                f"Message reference not found ({bug.message_snowflake})."
            )
        member = guild.get_member(member_snowflake)
        user_mentions = []
        for member_snowflake in bug.member_snowflakes:
            user_mentions.append(member.mention)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member.display_name} has been {action}",
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

    async def build_dictionary(self, where_kwargs):
        dictionary = {}
        bugs = await self.__database_factory.select(singular=False, **where_kwargs)
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

    async def build_pages(self, scope, where_kwargs, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Developer Logs"

        dictionary = await self.build_dictionary(where_kwargs=where_kwargs)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=BugDictionary, dictionary=dictionary
        )
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_messages)

        bug_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            field_count = 0
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for message_snowflake, entry in guild_data.get("messages", {}).items():
                channel = guild.get_channel(entry["channel_snowflake"])
                if not channel:
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
                if where_kwargs.get("id", None) == str(entry["id"]):
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
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name=f"**Channel:** {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                    lines = []
            if lines:
                embed.add_field(
                    name=f"**Channel:** {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({bug_n})**"
        return pages

    async def create_bug(self, message, reference):
        try:
            bug = Bug(
                channel_snowflake=message.channel.id,
                member_snowflakes=[],
                guild_snowflake=message.guild.id,
                id=reference,
                message_snowflake=message.id,
            )
            await self.__database_factory.create(bug)
        except discord.Forbidden as e:
            self.__bot.logger.info(str(e).capitalize())

    async def assign_bug_to_developer(
        self, developer, member_dict, reference, where_kwargs
    ):
        bug = await self.__database_factory.select(
            id=reference, resolved=False, singular=True
        )
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
