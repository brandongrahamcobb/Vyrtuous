"""!/bin/python3
streaming_service.py The purpose of this program is to extend Service service the stream command class.

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
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.stream import Stream
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.paginator import Paginator
from vyrtuous.utils.tracking.stream_embed import StreamEmbed
from vyrtuous.utils.users import moderator_service

MODEL = Stream


async def send_log(
    author_snowflake: int | None,
    channel_snowflake: int | None,
    duration_value: str | None,
    guild_snowflake: int,
    identifier: str,
    is_channel_scope: bool | None,
    member_snowflake: int,
    message_snowflake: int | None,
    message_channel_snowflake: int | None,
    role_snowflake: int | None,
    target: str | None,
    reason: str = "No reason provided",
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    embed = StreamEmbed(color=None, description=None, title=None, url=None)
    embed.set_title(identifier=identifier).set_action(
        duration_value=duration_value
    ).set_reason(reason=reason)
    pages: list[discord.Embed] = []
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    if channel_snowflake:
        if author_snowflake:
            executor_role = await moderator_service.resolve_highest_role(
                channel_snowflake=int(channel_snowflake),
                guild_snowflake=int(guild_snowflake),
                member_snowflake=int(author_snowflake),
            )
            embed.set_executor(
                author_snowflake=author_snowflake,
                guild_snowflake=guild_snowflake,
                highest_role=executor_role,
            )
            author = guild.get_member(author_snowflake)
            if author:
                embed.set_tn(url=author.display_avatar.url)
        else:
            executor_role = "Unknown"
        target_role = await moderator_service.resolve_highest_role(
            channel_snowflake=int(channel_snowflake),
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
        )
        embed.set_target(
            guild_snowflake=guild_snowflake,
            target=target,
            target_snowflake=member_snowflake,
            highest_role=target_role,
        )
        embed.set_description(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            target_snowflake=member_snowflake,
        )
        if is_channel_scope is not None:
            embed.set_channel_ctx(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                is_channel_scope=is_channel_scope,
            )
        if message_snowflake and message_channel_snowflake:
            await embed.set_message_ctx(
                guild_snowflake=guild_snowflake,
                identifier=identifier,
                message_snowflake=message_snowflake,
                message_channel_snowflake=message_channel_snowflake,
            )
            embed.set_reference(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                target_snowflake=member_snowflake,
                message_snowflake=message_snowflake,
            )
    if role_snowflake:
        embed.set_role(guild_snowflake=guild_snowflake, role_snowflake=role_snowflake)
    pages.append(embed)

    streaming = await database_factory.select(singular=False)
    for stream in streaming:
        channel_obj = bot.get_channel(stream.target_channel_snowflake)
        if channel_obj and isinstance(
            channel_obj,
            (discord.TextChannel, discord.VoiceChannel, discord.StageChannel),
        ):
            perms = channel_obj.permissions_for(channel_obj.guild.me)
            if perms.send_messages and not channel_obj.guild.me.is_timed_out():
                paginator = Paginator()
                await paginator.start_without_message(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=channel_obj.guild.id,
                    pages=pages,
                )
    return


async def toggle_stream(
    target_channel_snowflake: int,
    guild_snowflake: int,
    *,
    source_channel_snowflake: int,
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    source_channel = guild.get_channel(source_channel_snowflake)
    if source_channel is None:
        raise commands.ChannelNotFound(str(source_channel_snowflake))
    target_channel = guild.get_channel(target_channel_snowflake)
    if target_channel is None:
        raise commands.ChannelNotFound(str(target_channel_snowflake))
    if source_channel_snowflake:
        stream = await database_factory.select(
            guild_snowflake=guild_snowflake,
            source_channel_snowflake=source_channel_snowflake,
            target_channel_snowflake=target_channel_snowflake,
            singular=True,
        )
        if stream:
            await database_factory.delete(
                guild_snowflake=guild_snowflake,
                source_channel_snowflake=source_channel_snowflake,
                target_channel_snowflake=target_channel_snowflake,
            )
        source_text = f"from {source_channel_snowflake}"
    else:
        stream = await database_factory.select(
            guild_snowflake=guild_snowflake,
            target_channel_snowflake=target_channel_snowflake,
            singular=True,
        )
        if stream:
            await database_factory.delete(
                guild_snowflake=target_channel_snowflake,
                source_channel_snowflake=source_channel_snowflake,
            )
        source_text = "from all channels"
    if stream:
        action = "disabled"
    else:
        action = "enabled"
    embed = discord.Embed(
        title=f"{emojis.get_random_emoji()} Tracking {action.capitalize()} {source_text} to {target_channel.mention}.",
        color=0x00FF00,
    )
    return embed
