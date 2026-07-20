"""!/bin/python3
text_mute_service.py The purpose of this program is to extend AliasService to service the text mute infraction.

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

from datetime import datetime, timezone

import discord

from vyrtuous.aliases import text_mute_alias_service, untext_mute_alias_service
from vyrtuous.aliases.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.text_mute import TextMute

MODEL = TextMute


async def enforce_or_undo(
    alias_ctx: AliasContext, message: discord.Message
) -> discord.Embed:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    text_mute = await database_factory.select(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        member_snowflake=alias_ctx.member_snowflake,
        singular=True,
    )
    if text_mute:
        untext_mute_ctx = untext_mute_alias_service.UntextMuteMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
        )
        embed = await untext_mute_alias_service.untext_mute_by_message(
            ctx=untext_mute_ctx, display=True
        )
        return embed
    else:
        text_mute_ctx = text_mute_alias_service.TextMuteMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            duration=alias_ctx.duration,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            reason=alias_ctx.reason,
        )
        embed = await text_mute_alias_service.text_mute_by_message(
            ctx=text_mute_ctx, display=True
        )
        return embed


async def is_text_muted(
    channel: discord.abc.GuildChannel, member: discord.Member
) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    mute = await database_factory.select(
        channel_snowflake=channel.id,
        guild_snowflake=channel.guild.id,
        member_snowflake=member.id,
        singular=True,
    )
    if mute:
        return True
    return False


async def is_text_muted_then_mute_and_reset_cooldown(
    channel: discord.abc.GuildChannel, member: discord.Member
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    if await is_text_muted(channel=channel, member=member):
        try:
            await channel.set_permissions(
                member,
                send_messages=False,
                add_reactions=False,
                reason="Reinstating active text-mute.",
            )
            await update_last_text_muted(channel=channel, member=member)
        except discord.Forbidden as e:
            bot.logger.warning(e)


async def update_last_text_muted(
    channel: discord.abc.GuildChannel, member: discord.Member
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    where_kwargs = {
        "channel_snowflake": channel.id,
        "guild_snowflake": channel.guild.id,
        "member_snowflake": member.id,
    }
    set_kwargs = {
        "last_muted": datetime.now(timezone.utc),
        "reset": False,
    }
    await database_factory.update(
        set_kwargs=set_kwargs,
        where_kwargs=where_kwargs,
    )
    bot.logger.info(
        f"Updated last_muted record for {member.display_name} in {channel.name}."
    )
