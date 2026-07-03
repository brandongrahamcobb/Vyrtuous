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

from datetime import datetime, timedelta, timezone

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
    database_factory = DatabaseFactory(MODEL)
    text_mute = await database_factory.select(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        member_snowflake=alias_ctx.member_snowflake,
        singular=True,
    )
    if text_mute:
        ctx = untext_mute_alias_service.UntextMuteMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
        )
        embed = await untext_mute_alias_service.untext_mute_by_message(
            ctx=ctx, display=True
        )
        return embed
    else:
        ctx = text_mute_alias_service.TextMuteMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            duration_value=alias_ctx.duration_value,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            reason=alias_ctx.reason,
        )
        embed = await text_mute_alias_service.text_mute_by_message(
            ctx=ctx, display=True
        )
        return embed


async def clean_expired():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    expired_text_mutes = await database_factory.select(expired=True, singular=False)
    if expired_text_mutes:
        for expired_text_mute in expired_text_mutes:
            channel_snowflake = int(expired_text_mute.channel_snowflake)
            guild_snowflake = int(expired_text_mute.guild_snowflake)
            member_snowflake = int(expired_text_mute.member_snowflake)
            kwargs = {
                "channel_snowflake": channel_snowflake,
                "guild_snowflake": guild_snowflake,
                "member_snowflake": member_snowflake,
            }
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate guild {guild_snowflake}, cleaning up expired text-mute."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired text-mute."
                )
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), cleaning up expired text-mute."
                )
                continue
            await database_factory.delete(**kwargs)
            try:
                await channel.set_permissions(
                    target=member,
                    send_messages=None,
                    add_reactions=None,
                    reason="Cleaning up expired text-mute",
                )
            except discord.Forbidden as e:
                bot.logger.error(str(e).capitalize())
            except discord.HTTPException as e:
                bot.logger.error(f"HTTP error removing expired text mute: {e}")


async def clean_overwrites():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    now = datetime.now(timezone.utc)
    text_mutes = await database_factory.select(singular=False)
    for text_mute in text_mutes:
        channel_snowflake = int(text_mute.channel_snowflake)
        guild_snowflake = int(text_mute.guild_snowflake)
        member_snowflake = int(text_mute.member_snowflake)
        where_kwargs = {
            "channel_snowflake": channel_snowflake,
            "guild_snowflake": guild_snowflake,
            "member_snowflake": member_snowflake,
        }
        set_kwargs = {"reset": True}
        if not text_mute.reset and text_mute.last_muted < now - timedelta(weeks=1):
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                bot.logger.info(
                    f"Unable to locate guild {guild_snowflake} for removing overwrite."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                bot.logger.info(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                )
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                bot.logger.info(
                    f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                )
                continue
            try:
                await channel.set_permissions(
                    target=member,
                    overwrite=None,
                    reason="Resetting text-mute overwrite",
                )
            except discord.Forbidden as e:
                bot.logger.error(str(e).capitalize())
            await database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )


async def migrate(kwargs):
    database_factory = DatabaseFactory(MODEL)
    await database_factory.update(**kwargs)


async def is_text_muted(channel: discord.abc.GuildChannel, member: discord.Member):
    database_factory = DatabaseFactory(MODEL)
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
):
    bot = DiscordBot.get_instance()
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
):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
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
