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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.stream import stream_service
from vyrtuous.vegan.vegan import Vegan

MODEL = Vegan


async def enforce_or_undo(
    ctx,
    source: Union[commands.Context, discord.Interaction, discord.Message],
    tick,
):
    database_factory = DatabaseFactory(MODEL)
    obj = await database_factory.select(
        channel_snowflake=ctx.channel.id,
        guild_snowflake=ctx.guild.id,
        member_snowflake=ctx.member.id,
        singular=True,
    )
    if obj:
        await undo(ctx=ctx, source=source, tick=tick)
    else:
        await enforce(ctx=ctx, source=source, tick=tick)


async def enforce(ctx, source, tick):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    vegan = MODEL(
        guild_snowflake=ctx.guild.id,
        member_snowflake=ctx.member_snowflake,
    )
    await database_factory.create(vegan)
    member = ctx.guild.get_member(ctx.member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(
            ctx.member_snowflake, None
        )
        if simplified_member is None:
            raise commands.MemberNotFound(str(ctx.member_snowflake))
    await stream_service.send_log(
        channel=ctx.channel,
        identifier="vegan",
        member=member,
        source=source,
    )
    embed = await act_embed(ctx=ctx)
    return await tick.end(success=embed)


async def undo(ctx, source, tick):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    member = ctx.guild.get_member(ctx.member_snowflake)
    if member is None:
        simplified_member = bot.registry.get(MemberState).active.get(
            ctx.member_snowflake, None
        )
        if simplified_member is None:
            raise commands.MemberNotFound(str(ctx.member_snowflake))

    await database_factory.delete(
        channel_snowflake=ctx.channel.id,
        guild_snowflake=ctx.guild.id,
        member_snowflake=ctx.member_snowflake,
    )
    await stream_service.send_log(
        channel=ctx.channel,
        identifier="carnist",
        is_modification=True,
        member=member,
        source=source,
    )
    embed = await undo_embed(ctx=ctx)
    return await tick.end(success=embed)


async def act_embed(ctx):
    bot = DiscordBot.get_instance()
    member = ctx.guild.get_member(ctx.member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            ctx.member_snowflake, None
        )
        if simplified_member:
            display_name = simplified_member[0]
            member_str = display_name
        else:
            raise commands.MemberNotFound(str(ctx.member_snowflake))
    embed = discord.Embed(
        title=f"\U0001f525\U0001f525 {display_name} "
        f"is going Vegan!!!\U0001f525\U0001f525",
        description=(f"**User:** {member_str}\n"),
        color=discord.Color.blue(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed


async def undo_embed(ctx):
    bot = DiscordBot.get_instance()
    member = ctx.guild.get_member(ctx.member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            ctx.member_snowflake, None
        )
        if simplified_member:
            display_name = simplified_member[0]
            member_str = display_name
        else:
            raise commands.MemberNotFound(str(ctx.member_snowflake))
    embed = discord.Embed(
        title=f"\U0001f44e\U0001f44e "
        f"{display_name} is a Carnist \U0001f44e\U0001f44e",
        description=(f"**User:** {member_str}\n"),
        color=discord.Color.yellow(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed


async def migrate(kwargs):
    database_factory = DatabaseFactory(MODEL)
    await database_factory.update(**kwargs)
