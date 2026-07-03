"""!/bin/python3
ban_service.py The purpose of this program is to extend AliasService to service ban infractions.

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
from discord.ext import commands

from vyrtuous.aliases import ban_alias_service, unban_alias_service
from vyrtuous.aliases.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.ban import Ban
from vyrtuous.db.database_factory import DatabaseFactory

MODEL = Ban


async def enforce_or_undo(
    alias_ctx: AliasContext,
    message: discord.Message,
):
    database_factory = DatabaseFactory(MODEL)
    ban = await database_factory.select(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        member_snowflake=alias_ctx.member_snowflake,
        singular=True,
    )
    if ban:
        ctx = unban_alias_service.UnbanMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_channel_snowflake=message.channel.id,
            message_snowflake=message.id,
        )
        embed = await unban_alias_service.unban_by_message(ctx=ctx, display=True)
        return embed
    else:
        ctx = ban_alias_service.BanMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            duration_value=alias_ctx.duration_value,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            reason=alias_ctx.reason,
        )
        embed = await ban_alias_service.ban_by_message(
            ctx=ctx,
            display=True,
        )
        return embed


async def clean_expired():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    expired_bans = await database_factory.select(expired=True, singular=False)
    if expired_bans:
        for expired_ban in expired_bans:
            channel_snowflake = int(expired_ban.channel_snowflake)
            guild_snowflake = int(expired_ban.guild_snowflake)
            member_snowflake = int(expired_ban.member_snowflake)
            kwargs = {
                "channel_snowflake": channel_snowflake,
                "guild_snowflake": guild_snowflake,
                "member_snowflake": member_snowflake,
            }
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate guild {guild_snowflake}, cleaning up expired ban."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}, cleaning up expired ban."
                )
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), cleaning up expired ban."
                )
                continue
            await database_factory.delete(**kwargs)
            try:
                await channel.set_permissions(
                    member, view_channel=None, reason="Cleaning up expired ban."
                )
            except discord.Forbidden as e:
                bot.logger.error(str(e).capitalize())
            except discord.HTTPException as e:
                bot.logger.error(f"HTTP error removing expired ban: {e}")


async def toggle_blacklist(channel, member_snowflake: int):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    ban = await database_factory.select(
        channel_snowflake=channel.id,
        member_snowflake=member_snowflake,
        singular=True,
    )
    member = channel.guild.get_member(member_snowflake)
    if member:
        display_name = member.display_name
    else:
        member_data = bot.registry.get(MemberState).active.get(member_snowflake, None)
        if member_data:
            display_name = member_data[0]
        else:
            raise commands.MemberNotFound(str(member_snowflake))
    if not ban:
        return f"{display_name} is not banned in {channel.mention}."
    where_kwargs = {
        "channel_snowflake": channel.id,
        "member_snowflake": member_snowflake,
    }
    if ban.blacklisted:
        set_kwargs = {"blacklisted": False}
        action = "unlisted"
        await database_factory.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
    else:
        set_kwargs = {"blacklisted": True}
        action = "blacklisted"
        await database_factory.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
    bot.logger.info(
        f"{display_name} ({member_snowflake}) has been ban {action} in {channel.mention}."
    )
    return f"{display_name} ({member_snowflake}) has been ban {action} in {channel.mention}."


async def clean_overwrites():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    bans = await database_factory.select(singular=False)
    for ban in bans:
        channel_snowflake = int(ban.channel_snowflake)
        guild_snowflake = int(ban.guild_snowflake)
        member_snowflake = int(ban.member_snowflake)
        where_kwargs = {
            "channel_snowflake": channel_snowflake,
            "guild_snowflake": guild_snowflake,
            "member_snowflake": member_snowflake,
        }
        set_kwargs = {"reset": True}
        if (
            not ban.reset
            and ban.last_kicked < datetime.now(timezone.utc) - timedelta(weeks=1)
            and not ban.blacklisted
        ):
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
                    target=member, overwrite=None, reason="Resetting ban overwrite."
                )
            except discord.Forbidden as e:
                bot.logger.error(str(e).capitalize())
            await database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )


async def migrate(kwargs):
    database_factory = DatabaseFactory(MODEL)
    await database_factory.update(**kwargs)


async def is_banned(channel: discord.abc.GuildChannel, member: discord.Member):
    database_factory = DatabaseFactory(MODEL)
    ban = await database_factory.select(
        channel_snowflake=channel.id, member_snowflake=member.id, singular=True
    )
    if ban:
        return True
    return False


async def is_banned_then_kick_and_reset_cooldown(
    channel: discord.abc.GuildChannel, member: discord.Member
):
    if await is_banned(channel=channel, member=member):
        if (
            member.voice
            and member.voice.channel
            and member.voice.channel.id == channel.id
        ):
            await kick(channel=channel, member=member)
        await toggle_view_channel(channel=channel, member=member, view_channel=False)


async def toggle_view_channel(
    channel: discord.abc.GuildChannel,
    member: discord.Member,
    view_channel: bool,
):
    bot = DiscordBot.get_instance()
    try:
        await channel.set_permissions(
            member,
            view_channel=view_channel,
            reason=f"Toggled ban {'off' if not view_channel else 'on'}.",
        )
    except discord.Forbidden as e:
        bot.logger.warning(e)


async def kick(channel: discord.abc.GuildChannel, member: discord.Member):
    bot = DiscordBot.get_instance()
    try:
        await member.move_to(None, reason="Reinstating active ban.")
        await update_last_kicked(channel=channel, member=member)
    except discord.Forbidden as e:
        bot.logger.warning(e)


async def update_last_kicked(channel: discord.abc.GuildChannel, member: discord.Member):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    where_kwargs = {"channel_snowflake": channel.id, "member_snowflake": member.id}
    set_kwargs = {
        "last_kicked": datetime.now(timezone.utc),
        "reset": False,
    }
    await database_factory.update(
        set_kwargs=set_kwargs,
        where_kwargs=where_kwargs,
    )
    bot.logger.info(
        f"Updated last_kicked record for {member.display_name} in {channel.name}."
    )
