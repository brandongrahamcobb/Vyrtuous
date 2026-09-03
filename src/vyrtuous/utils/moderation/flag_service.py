"""!/bin/python3
flag_service.py The purpose of this program is to service flags.

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
) -> discord.Embed | str:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    obj = await database_factory.select(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        member_snowflake=alias_ctx.member_snowflake,
        singular=True,
    )
    if obj:
        unflag_ctx = unflag_alias_service.UnflagMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
        )
        embed = await unflag_alias_service.unflag_by_message(
            ctx=unflag_ctx, display=True
        )
        return embed
    else:
        flag_ctx = flag_alias_service.FlagMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            reason=alias_ctx.reason,
        )
        embed_or_message = await flag_alias_service.flag_by_message(
            ctx=flag_ctx, display=True
        )
        return embed_or_message


async def warn(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        return
    member = guild.get_member(member_snowflake)
    if member is None:
        return
    flags = bot.registry.get(MemberState).flagged
    flag_reason = flags[guild_snowflake].get(member.id, {}).get(channel.id, None)
    if flag_reason is not None and "Vegan" in channel.name:
        if isinstance(channel, discord.channel.VocalGuildChannel):
            if bot.registry.get(ChannelState).should_notify(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                member_snowflake=member_snowflake,
                timeout=300.0,
            ):
                embed = discord.Embed(
                    title=f"\u26a0\ufe0f {member.display_name} is flagged",
                    description=f"**Channel:** {channel.mention}\n**Reason:** {flag_reason}",
                    color=discord.Color.red(),
                )
                embed.set_thumbnail(url=member.display_avatar.url)
                await channel.send(embed=embed)
        bot.registry.get(ChannelState).record(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
        )


async def populate() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    original_dict = bot.registry.get(MemberState).flagged
    flags = await database_factory.select(singular=False)
    for flag in flags:
        original_dict[flag.guild_snowflake].setdefault(
            flag.member_snowflake,
            {}
        )[flag.channel_snowflake] = flag.reason
