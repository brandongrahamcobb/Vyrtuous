"""!/bin/python3
flag_service.py The purpose of this program is to extend AliasService to service flag infractions.

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

import time
from copy import copy
from dataclasses import dataclass, field
from typing import Dict, List, Union

import discord
from discord.ext import commands

from vyrtuous.flag.flag import Flag


@dataclass
class FlagDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, str]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class FlagService:
    __CHUNK_SIZE = 7
    MODEL = Flag

    def __init__(
        self,
        *,
        bot,
        database_factory=None,
        data_service=None,
        dictionary_service=None,
        emoji=None,
        stream_service=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__data_service = data_service
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__flags = []
        self.__join_log = {}
        self.__stream_service = stream_service

    async def enforce_or_undo(
        self,
        ctx,
        source: Union[commands.Context, discord.Interaction, discord.Message],
        state,
    ):
        obj = await self.__database_factory.select(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            singular=True,
        )
        if obj:
            await self.undo(ctx=ctx, source=source, state=state)
        else:
            await self.enforce(ctx=ctx, source=source, state=state)

    async def build_dictionary(self, where_kwargs):
        dictionary = {}
        flags = await self.__database_factory.select(singular=False, **where_kwargs)
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
        return dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        bot = self.__bot
        title = f"{self.__emoji.get_random_emoji()} Flags {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)

        dictionary = await self.build_dictionary(where_kwargs=where_kwargs)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=FlagDictionary, dictionary=dictionary
        )
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)

        flag_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, flag_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
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
                        lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
                    flag_n += 1
                    field_count += 1
                    if field_count >= self.__CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(lines),
                            inline=False,
                        )
                        embed = self.__dictionary_service.flush_page(
                            embed, pages, title, guild.name
                        )
                        lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({flag_n})**"
        return pages

    async def enforce(self, ctx, source, state):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        author = guild.get_member(ctx.author_snowflake)
        channel = guild.get_channel(ctx.target_channel_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        flag = self.MODEL(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            reason=ctx.reason,
        )
        await self.__database_factory.create(flag)
        self.__flags.append(flag)
        await self.__stream_service.send_log(
            channel=channel,
            identifier="flag",
            member=member,
            source=source,
            reason=ctx.reason,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="flag",
            reason=ctx.reason,
            target=member,
        )
        embed = await self.act_embed(ctx=ctx)
        return await state.end(success=embed)

    async def undo(self, ctx, source, state):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        author = guild.get_member(ctx.author_snowflake)
        channel = guild.get_channel(ctx.target_channel_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        await self.__database_factory.delete(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
        )
        for flag in self.__flags:
            if flag.channel_snowflake == ctx.target_channel_snowflake:
                self.__flags.remove(flag)
                break
        await self.__stream_service.send_log(
            channel=channel,
            identifier="unflag",
            is_modification=True,
            member=member,
            source=source,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="uflag",
            is_modification=True,
            reason=ctx.reason,
            target=member,
        )
        embed = await self.undo_embed(ctx=ctx)
        return await state.end(success=embed)

    async def act_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        channel = guild.get_channel(ctx.target_channel_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member.display_name} has been flagged",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Reason:** {ctx.reason}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def undo_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        channel = guild.get_channel(ctx.target_channel_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member.display_name} has been unflagged",
            description=(f"**User:** {member.mention}\n**Channel:** {channel.mention}"),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def migrate(self, kwargs):
        self.__database_factory.update(**kwargs)

    async def load_flags_into_memory(self):
        self.__flags = await self.__database_factory.select()

    async def warn(self, channel: discord.abc.GuildChannel, member: discord.Member):
        if channel.id == 1222056499959042108:
            for flag in self.__flags:
                if flag.channel_snowflake == channel.id:
                    if flag.member_snowflake == member.id:
                        embed = discord.Embed(
                            title=f"\u26a0\ufe0f {member.display_name} is flagged",
                            color=discord.Color.red(),
                        )
                        embed.set_thumbnail(url=member.display_avatar.url)
                        embed.add_field(
                            name=f"Channel: {channel.mention}",
                            value=f"Reason: {flag.reason}",
                            inline=False,
                        )
                        now = time.time()
                        self.__join_log[member.id] = [
                            t for t in self.__join_log[member.id] if now - t < 300
                        ]
                        if len(self.__join_log[member.id]) < 1:
                            self.__join_log[member.id].append(now)
                            await channel.send(embed=embed)
