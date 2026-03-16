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
    __CHUNK_SIZE = 12
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
        **kwargs,
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
        default_ctx,
        source: Union[commands.Context, discord.Interaction, discord.Message],
        state,
    ):
        obj = await self.__database_factory.select(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member.id,
            singular=True,
        )
        if obj:
            await self.undo(
                ctx=ctx, default_ctx=default_ctx, source=source, state=state
            )
        else:
            await self.enforce(
                ctx=ctx, default_ctx=default_ctx, source=source, state=state
            )

    async def build_dictionary(self, obj):
        flags = []
        dictionary = {}
        if isinstance(obj, discord.Guild):
            flags = await self.__database_factory.select(guild_snowflake=obj.id)
        elif isinstance(obj, discord.abc.GuildChannel):
            flags = await self.__database_factory.select(channel_snowflake=obj.id)
        elif isinstance(obj, discord.Member):
            flags = await self.__database_factory.select(member_snowflake=obj.id)
        else:
            flags = await self.__database_factory.select()
        if flags:
            for flag in flags:
                dictionary.setdefault(flag.guild_snowflake, {"members": {}})
                dictionary[flag.guild_snowflake]["members"].setdefault(
                    flag.member_snowflake, {"flags": {}}
                )
                dictionary[flag.guild_snowflake]["members"][flag.member_snowflake][
                    "flags"
                ].setdefault(flag.channel_snowflake, {})
                dictionary[flag.guild_snowflake]["members"][flag.member_snowflake][
                    "flags"
                ][flag.channel_snowflake].update(
                    {
                        "reason": flag.reason,
                    }
                )
        return dictionary

    async def build_pages(self, is_at_home, obj):
        lines, pages = [], []

        obj_name = "All Servers"
        if obj:
            obj_name = obj.name
        title = f"{self.__emoji.get_random_emoji()} Flags for {obj_name}"

        dictionary = await self.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=FlagDictionary, dictionary=dictionary
        )

        for guild_snowflake, guild_data in processed_dictionary.data.items():
            flag_n = 0
            field_count = 0
            lines = []
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, flag_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(obj, discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(url=obj.display_avatar.url)
                    thumbnail = True
                for channel_snowflake, channel_dictionary in flag_dictionary.get(
                    "flags", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(obj, discord.abc.GuildChannel):
                        lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(obj, discord.Member):
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
            original_description = embed.description or ""
            embed.description = f"**{original_description} ({flag_n})**"
            pages.append(embed)
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)
        return pages

    async def delete(
        self,
        author,
        source,
        *,
        guild_snowflake=None,
        channel_snowflake=None,
        member_snowflake=None,
    ):
        kwargs = {}
        if channel_snowflake:
            kwargs.update({"channel_snowflake": channel_snowflake})
        if guild_snowflake:
            kwargs.update({"guild_snowflake": guild_snowflake})
        if member_snowflake:
            kwargs.update({"member_snowflake": member_snowflake})
        objects = await self.__database_factory.select(**kwargs)
        for obj in objects:
            await self.__database_factory.delete_by_cls(obj, **kwargs)
            guild = self.__bot.get_guild(obj.guild_snowflake)
            channel = guild.get_channel(obj.channel_snowflake)
            member = guild.get_member(obj.member_snowflake)
            await self.undo_log(
                author=author,
                channel=channel,
                member=member,
                source=source,
            )

    async def enforce_log(self, author, channel, member, reason, source):
        await self.__stream_service.send_log(
            author=author,
            channel=channel,
            identifier="flag",
            member=member,
            reason=reason,
            source=source,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="flag",
            member=member,
            reason=reason,
        )

    async def enforce(self, ctx, default_ctx, source, state):
        flag = self.MODEL(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member.id,
            reason=ctx.reason,
        )
        await self.__database_factory.create(flag)
        self.__flags.append(flag)
        await self.enforce_log(
            author=default_ctx.author,
            channel=ctx.channel,
            member=ctx.member,
            reason=ctx.reason,
            source=source,
        )
        embed = await self.act_embed(ctx=ctx)
        return await state.end(success=embed)

    async def undo_log(self, author, channel, member, source):
        await self.__stream_service.send_log(
            author=author,
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
            member=member,
        )

    async def undo(self, ctx, default_ctx, source, state):
        await self.__database_factory.delete(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member.id,
        )
        for flag in self.__flags:
            if flag.channel_snowflake == ctx.channel.id:
                self.__flags.remove(flag)
                break
        await self.undo_log(
            author=default_ctx.author,
            channel=ctx.channel,
            member=ctx.member,
            source=source,
        )
        embed = await self.undo_embed(ctx=ctx)
        return await state.end(success=embed)

    async def act_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.guild.id)
        channel = guild.get_channel(ctx.channel.id)
        member = guild.get_member(ctx.member.id)
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
        guild = self.__bot.get_guild(ctx.guild.id)
        channel = guild.get_channel(ctx.channel.id)
        member = guild.get_member(ctx.member.id)
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
                            name=f"User: {member.mention}",
                            value=f"Channel: {channel.mention}\nReason: {flag.reason}",
                            inline=False,
                        )
                        now = time.time()
                        self.__join_log[member.id] = [
                            t for t in self.__join_log[member.id] if now - t < 300
                        ]
                        if len(self.__join_log[member.id]) < 1:
                            self.__join_log[member.id].append(now)
                            await channel.send(embed=embed)
