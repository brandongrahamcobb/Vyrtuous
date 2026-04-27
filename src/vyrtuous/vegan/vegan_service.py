"""!/bin/python3
vegan_service.py The purpose of this program is to extend Service to service the vegan class.

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
from typing import Dict, List, Union

from discord.ext import commands
import discord


from vyrtuous.active_members import active_member_service
from vyrtuous.vegan.vegan import Vegan


@dataclass
class VeganDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, dict]]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class VeganService:
    __CHUNK_SIZE = 12
    MODEL = Vegan
    vegans = {}

    def __init__(
        self,
        *,
        active_member_service=None,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
        stream_service=None,
    ):
        self.__active_member_service = active_member_service
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__stream_service = stream_service

    async def populate(self):
        vegans = await self.__database_factory.select()
        for vegan in vegans:
            guild = self.__bot.get_guild(vegan.guild_snowflake)
            if not guild:
                continue
            self.vegans[vegan.member_snowflake] = {
                "last_active": None,
                "name": vegan.display_name,
            }

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

    async def build_clean_dictionary(self, obj):
        vegans = []
        dictionary = {}
        if isinstance(obj, discord.Guild):
            vegans = await self.__database_factory.select(guild_snowflake=obj.id)
        elif isinstance(obj, discord.abc.GuildChannel):
            vegans = await self.__database_factory.select(channel_snowflake=obj.id)
        elif isinstance(obj, discord.Member):
            vegans = await self.__database_factory.select(member_snowflake=obj.id)
        else:
            vegans = await self.__database_factory.select()
        if vegans:
            for vegan in vegans:
                dictionary.setdefault(vegan.guild_snowflake, {"members": {}})
                dictionary[vegan.guild_snowflake]["members"].setdefault(
                    vegan.member_snowflake, {"vegans": {}}
                )
                dictionary[vegan.guild_snowflake]["members"][vegan.member_snowflake][
                    "vegans"
                ].setdefault("placeholder", {})
        return dictionary

    async def build_pages(self, is_at_home, obj):
        lines, pages = [], []

        obj_name = "All Servers"
        if not isinstance(obj, int):
            obj_name = obj.name
        else:
            member = self.__active_member_service.active_member.get(obj, None)
            if member:
                obj_name = member.get("name", None)
            else:
                return "No vegans found."
        title = f"{self.__emoji.get_random_emoji()} Vegans for {obj_name}"

        dictionary = await self.__dictionary_service.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=VeganDictionary, dictionary=dictionary
        )

        vegan_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            field_count = 0
            lines = []
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, vegan_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not member:
                    if not isinstance(obj, discord.Member):
                        lines.append(
                            f"**User:** {member.display_name} {member.mention}"
                        )
                    else:
                        if not thumbnail:
                            embed.set_thumbnail(url=obj.display_avatar.url)
                            thumbnail = True
                else:
                    member = self.__active_member_service.active_members.get(
                        member_snowflake, None
                    )
                    if member:
                        display_name = member.get("name", None)
                        lines.append(f"**User:** {display_name} ({member_snowflake})")
                vegan_n += 1
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
                    field_count = 0
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)
        if pages:
            original_description = embed.description or ""
            embed.description = f"**{original_description} ({vegan_n})**"
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)
        return pages

    async def enforce(self, ctx, source, state):
        vegan = self.MODEL(
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member_snowflake,
        )
        await self.__database_factory.create(vegan)
        member = ctx.guild.get_member(ctx.member_snowflake)
        if not member:
            member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
        self.vegans.update({ctx.member_snowflake: {"name": ctx.display_name}})
        await self.__stream_service.send_entry(
            channel_snowflake=ctx.channel.id,
            identifier="vegan",
            member=member,
            source=source,
        )
        embed = await self.act_embed(ctx=ctx)
        return await state.end(success=embed)

    async def undo(self, ctx, source, state):
        member = ctx.guild.get_member(ctx.member.id)
        if not member:
            member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
        del self.vegans[ctx.member_snowflake]
        await self.__database_factory.delete(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member_snowflake,
        )
        await self.__stream_service.send_entry(
            channel_snowflake=ctx.channel.id,
            identifier="carnist",
            is_modification=True,
            member=member,
            source=source,
        )
        embed = await self.undo_embed(ctx=ctx)
        return await state.end(success=embed)

    async def act_embed(self, ctx):
        member = ctx.guild.get_member(ctx.member.id)
        if member:
            member_display_name = member.display_name
            member_str = member.mention
        else:
            simplified_member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
            member_display_name = simplified_member.get("name", None)
            member_str = simplified_member.get("name", None)
        embed = discord.Embed(
            title=f"\U0001f525\U0001f525 {member_display_name} "
            f"is going Vegan!!!\U0001f525\U0001f525",
            description=(f"**User:** {member_str}\n"),
            color=discord.Color.blue(),
        )
        if member:
            embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def undo_embed(self, ctx):
        member = ctx.guild.get_member(ctx.member.id)
        if member:
            member_display_name = member.display_name
            member_str = member.mention
        else:
            simplified_member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
            member_display_name = simplified_member.get("name", None)
            member_str = simplified_member.get("name", None)
        embed = discord.Embed(
            title=f"\U0001f44e\U0001f44e "
            f"{member_display_name} is a Carnist \U0001f44e\U0001f44e",
            description=(f"**User:** {member_str}\n"),
            color=discord.Color.yellow(),
        )
        if member:
            embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def migrate(self, kwargs):
        self.__database_factory.update(**kwargs)
