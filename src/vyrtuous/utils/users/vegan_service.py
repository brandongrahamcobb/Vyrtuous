"""!/bin/python3
vegan_service.py The purpose of this program is to service the vegans.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState, MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.vegan import Vegan
from vyrtuous.utils.messaging import emojis

MODEL = Vegan


def is_vegan(guild_snowflake: int, member_snowflake: int) -> bool:
    bot: DiscordBot = DiscordBot.get_instance()
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


async def toggle_vegan(
    guild_snowflake: int, member_snowflake: int, notes: str | None
) -> discord.Embed:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
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
        embed = await build_carnist_embed(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )
        return embed
    else:
        vegan = MODEL(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            notes=notes,
        )
        await database_factory.create(vegan)
        embed = await build_vegan_embed(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            notes=notes,
        )
        return embed


async def build_vegan_embed(
    guild_snowflake: int, member_snowflake: int, notes: str | None
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
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
    bot: DiscordBot = DiscordBot.get_instance()
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


async def notify(
    channel: discord.channel.VocalGuildChannel, member: discord.Member
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    vegans = await database_factory.select(singular=False)
    for vegan in vegans:
        if "Vegan" in channel.name and vegan.member_snowflake == member.id:
            if bot.registry.get(ChannelState).should_notify(channel.id, member.id):
                embed = discord.Embed(
                    title=f"{emojis.get_random_emoji()}) {member.display_name} is a recent Vegan!",
                    description=f"**Notes:** {vegan.notes}",
                    color=discord.Color.green(),
                )
                embed.set_thumbnail(url=member.display_avatar.url)
                await channel.send(embed=embed)
            bot.registry.get(ChannelState).record(channel.id, member.id)
