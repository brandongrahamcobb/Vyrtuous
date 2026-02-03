"""vegan.py The purpose of this program is to inherit from the DatabaseFactory to provide the vegan role.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.service.mgmt.alias_service import AliasService
from vyrtuous.service.mgmt.stream_service import StreamService
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji


class VeganService(AliasService):

    lines, pages = [], []
    model = Vegan

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        vegans = await Vegan.select(singular=False, **where_kwargs)
        for vegan in vegans:
            dictionary.setdefault(vegan.guild_snowflake, {"members": {}})
            dictionary[vegan.guild_snowflake]["members"].setdefault(
                vegan.member_snowflake, {"vegans": {}}
            )
            dictionary[vegan.guild_snowflake]["members"][vegan.member_snowflake][
                "vegans"
            ].setdefault("placeholder", {})
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                VeganService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                VeganService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Vegans {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await VeganService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, vegan_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    VeganService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                else:
                    if not thumbnail:
                        embed.set_thumbnail(
                            url=object_dict.get("object", None).display_avatar.url
                        )
                        thumbnail = True
                field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(VeganService.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, VeganService.pages, title, guild.name)
                    VeganService.lines = []
                    field_count = 0
            if VeganService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(VeganService.lines),
                    inline=False,
                )
            VeganService.pages.append(embed)
        return VeganService.pages

    @classmethod
    async def enforce(cls, information, message, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        vegan = Vegan(
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
        )
        await vegan.create()
        await StreamService.send_entry(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            identifier="vegan",
            member=member,
            message=message,
        )
        embed = await VeganService.act_embed(information=information)
        return await state.end(success=embed)

    @classmethod
    async def undo(cls, information, message, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        await Vegan.delete(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
        )
        await StreamService.send_entry(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            identifier="carnist",
            is_modification=True,
            member=member,
            message=message,
        )
        embed = await VeganService.undo_embed(information=information)
        return await state.end(success=embed)

    @classmethod
    async def act_embed(cls, information):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"\U0001f525\U0001f525 {member.display_name} "
            f"is going Vegan!!!\U0001f525\U0001f525",
            description=(f"**User:** {member.mention}\n"),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, information):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"\U0001f44e\U0001f44e "
            f"{member.display_name} is a Carnist \U0001f44e\U0001f44e",
            description=(f"**User:** {member.mention}\n"),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
