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

import discord

from vyrtuous.aliases import flag_alias_service, unflag_alias_service
from vyrtuous.aliases.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState, MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.flag import Flag

MODEL = Flag


async def enforce_or_undo(
    alias_ctx: AliasContext,
    message: discord.Message,
) -> discord.Embed:
    database_factory = DatabaseFactory(MODEL)
    obj = await database_factory.select(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        member_snowflake=alias_ctx.member_snowflake,
        singular=True,
    )
    if obj:
        ctx = flag_alias_service.FlagMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            reason=alias_ctx.reason,
        )
        embed = await flag_alias_service.flag_by_message(ctx=ctx, display=True)
        return embed
    else:
        ctx = unflag_alias_service.UnflagMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
        )
        embed = await unflag_alias_service.unflag_by_message(ctx=ctx, display=True)
        return embed


async def migrate(kwargs):
    database_factory = DatabaseFactory(MODEL)
    await database_factory.update(**kwargs)


async def warn(channel: discord.channel.VocalGuildChannel, member: discord.Member):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    flags = await database_factory.select(singular=False)
    for flag in flags:
        if flag.channel_snowflake == channel.id and flag.member_snowflake == member.id:
            if bot.registry.get(ChannelState).should_warn(channel.id, member.id):
                embed = discord.Embed(
                    title=f"\u26a0\ufe0f {member.display_name} is flagged",
                    description=f"Channel: {channel.mention}\nReason: {flag.reason}",
                    color=discord.Color.red(),
                )
                embed.set_thumbnail(url=member.display_avatar.url)
                await channel.send(embed=embed)
            bot.registry.get(ChannelState).record(channel.id, member.id)


async def populate():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    original_set = bot.registry.get(MemberState).flagged
    flags = await database_factory.select(singular=False)
    for flag in flags:
        original_set[flag.guild_snowflake].add(flag.member_snowflake)
