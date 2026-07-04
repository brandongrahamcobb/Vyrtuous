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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.vegan import Vegan

MODEL = Vegan


def is_vegan(guild_snowflake: int, member_snowflake: int):
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    member = guild.get_member(member_snowflake)
    if member is None:
        raise commands.MemberNotFound(str(member_snowflake))
    for role in member.roles:
        if role.name == "Vegan":
            return True
    else:
        return False


async def toggle_vegan(guild_snowflake: int, member_snowflake: int, notes: str | None):
    database_factory = DatabaseFactory(MODEL)
    vegan = await database_factory.select(
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        singular=True,
    )
    if vegan:
        await database_factory.delete(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        embed = await build_vegan_embed(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            notes=notes,
        )
        return embed
    else:
        vegan = MODEL(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            notes=notes,
        )
        await database_factory.create(vegan)
        embed = await build_carnist_embed(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        return embed


async def build_vegan_embed(
    guild_snowflake: int, member_snowflake: int, notes: str | None
):
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member:
            display_name = simplified_member[0]
            member_str = display_name
        else:
            raise commands.MemberNotFound(str(member_snowflake))
    embed = discord.Embed(
        title=f"\U0001f525\U0001f525 {display_name} "
        f"is going Vegan!!!\U0001f525\U0001f525",
        description=(f"**User:** {member_str}\n**Notes:** {notes}\n"),
        color=discord.Color.blue(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed


async def build_carnist_embed(guild_snowflake: int, member_snowflake: int):
    bot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member:
            display_name = simplified_member[0]
            member_str = display_name
        else:
            raise commands.MemberNotFound(str(member_snowflake))
    embed = discord.Embed(
        title=f"\U0001f44e\U0001f44e "
        f"{display_name} is a Carnist \U0001f44e\U0001f44e",
        description=(f"**User:** {member_str}\n"),
        color=discord.Color.yellow(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed
