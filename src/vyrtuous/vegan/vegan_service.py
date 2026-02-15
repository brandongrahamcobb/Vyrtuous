from copy import copy

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

import discord

from vyrtuous.vegan.vegan import Vegan


class VeganService:
    __CHUNK_SIZE = 7
    MODEL = Vegan

    def __init__(
        self,
        *,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
        stream_service=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__stream_service = stream_service

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        pages = []
        dictionary = {}
        vegans = await self.__database_factory.select(singular=False, **where_kwargs)
        for vegan in vegans:
            dictionary.setdefault(vegan.guild_snowflake, {"members": {}})
            dictionary[vegan.guild_snowflake]["members"].setdefault(
                vegan.member_snowflake, {"vegans": {}}
            )
            dictionary[vegan.guild_snowflake]["members"][vegan.member_snowflake][
                "vegans"
            ].setdefault("placeholder", {})
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        skipped_members = self.__dictionary_service.generate_skipped_members(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        # if is_at_home:
        #     if skipped_guilds:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_set_pages(
        #                 skipped=skipped_guilds,
        #                 title="Skipped Servers",
        #             )
        #         )
        #     if skipped_members:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_dict_pages(
        #                 skipped=skipped_members,
        #                 title="Skipped Members in Server",
        #             )
        #         )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Vegans {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.__dictionary_service.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        vegan_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, vegan_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                else:
                    if not thumbnail:
                        embed.set_thumbnail(
                            url=object_dict.get("object", None).display_avatar.url
                        )
                        thumbnail = True
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
            pages[0].description = f"**({vegan_n})**"
        return pages

    async def enforce(self, ctx, source, state):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        vegan = self.MODEL(
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
        )
        await self.__database_factory.create(vegan)
        await self.__stream_service.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="vegan",
            member=member,
            source=source,
        )
        embed = await self.act_embed(ctx=ctx)
        return await state.end(success=embed)

    async def undo(self, ctx, source, state):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        await self.__database_factory.delete(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
        )
        await self.__stream_service.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="carnist",
            is_modification=True,
            member=member,
            source=source,
        )
        embed = await self.undo_embed(ctx=ctx)
        return await state.end(success=embed)

    async def act_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"\U0001f525\U0001f525 {member.display_name} "
            f"is going Vegan!!!\U0001f525\U0001f525",
            description=(f"**User:** {member.mention}\n"),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def undo_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"\U0001f44e\U0001f44e "
            f"{member.display_name} is a Carnist \U0001f44e\U0001f44e",
            description=(f"**User:** {member.mention}\n"),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
