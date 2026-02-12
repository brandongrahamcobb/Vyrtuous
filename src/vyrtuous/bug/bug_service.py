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

import discord

from vyrtuous.bug.bug import Bug


class BugService:
    __CHUNK_SIZE = 7
    MODEL = Bug

    def __init__(
        self, bot=None, database_factory=None, dictionary_service=None, emoji=None
    ):
        self.bot = bot
        self.dictionary_service = dictionary_service
        self.emoji = emoji
        self.database_factory = database_factory

    async def create_embed(self, action, bug, member_snowflake):
        guild = self.bot.get_guild(bug.guild_snowflake)
        channel = guild.get_channel(bug.channel_snowflake)
        try:
            msg = await channel.fetch_message(bug.message_snowflake)
        except discord.NotFound:
            self.bot.logger.warning(
                f"Message reference not found ({bug.message_snowflake})."
            )
        member = guild.get_member(member_snowflake)
        user_mentions = []
        for member_snowflake in bug.member_snowflakes:
            user_mentions.append(member.mention)
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member.display_name} has been {action}",
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

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        pages = []
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
        skipped_guilds = self.dictionary_service.generate_skipped_guilds(dictionary)
        skipped_messages = await self.dictionary_service.generate_skipped_messages(
            dictionary
        )
        cleaned_dictionary = self.dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_messages=skipped_messages,
        )
        if is_at_home:
            if skipped_guilds:
                pages.extend(
                    self.dictionary_service.generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_messages:
                pages.extend(
                    self.dictionary_service.generate_skipped_dict_pages(
                        skipped=skipped_messages,
                        title="Skipped Messages in Server",
                    )
                )
        return cleaned_dictionary

    async def build_pages(self, filter, where_kwargs, is_at_home):
        lines, pages = [], []
        title = f"{self.emoji.get_random_emoji()} Developer Logs"

        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        bug_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for message_snowflake, entry in guild_data.get("messages", {}).items():
                channel = guild.get_channel(entry["channel_snowflake"])
                if not channel:
                    continue
                if filter == "resolved" and not entry.get("resolved"):
                    continue
                if filter == "unresolved" and entry.get("resolved"):
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
                    embed = self.dictionary_service.flush_page(
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

    async def assign_bug_to_developer(self, reference, member_dict):
        where_kwargs = member_dict.get("columns", None)

        self.database_factory.model = Developer
        developer = await self.database_factory.select(singular=True, **where_kwargs)
        if not developer:
            return (
                f"Developer not found for target ({member_dict.get('mention', None)})."
            )
        self.database_factory.model = self.MODEL
        bug = await self.database_factory.select(
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
