"""!/bin/python3
voice_mute_service.py The purpose of this program is to extend AliasService to service the voice mute infraction.

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

from vyrtuous.aliases import unvoice_mute_alias_service, voice_mute_alias_service
from vyrtuous.aliases.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.users import moderator_service

MODEL = VoiceMute


async def enforce_or_undo(
    alias_ctx: AliasContext,
    message: discord.Message,
) -> discord.Embed:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    automute_database_factory: DatabaseFactory = DatabaseFactory(AutoMute)
    auto_mute_channel = await automute_database_factory.select(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        singular=True,
    )
    if auto_mute_channel:
        target = "auto"
        auto_mute = await database_factory.select(
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            target=target,
            singular=True,
        )
        if auto_mute:
            unvoice_mute_ctx = unvoice_mute_alias_service.UnvoiceMuteMessageContext(
                author_snowflake=message.author.id,
                channel_snowflake=alias_ctx.channel_snowflake,
                guild_snowflake=alias_ctx.guild_snowflake,
                member_snowflake=alias_ctx.member_snowflake,
                message_snowflake=message.id,
                message_channel_snowflake=message.channel.id,
                target=target,
            )
            embed = await unvoice_mute_alias_service.unvoice_mute_by_message(
                ctx=unvoice_mute_ctx, display=True
            )
            return embed
        else:
            voice_mute_ctx = voice_mute_alias_service.VoiceMuteMessageContext(
                author_snowflake=message.author.id,
                channel_snowflake=alias_ctx.channel_snowflake,
                duration_value=alias_ctx.duration_value,
                guild_snowflake=alias_ctx.guild_snowflake,
                member_snowflake=alias_ctx.member_snowflake,
                message_snowflake=message.id,
                message_channel_snowflake=message.channel.id,
                reason=alias_ctx.reason,
                target=target,
            )
            embed = await voice_mute_alias_service.voice_mute_by_message(
                ctx=voice_mute_ctx, display=True
            )
            return embed
    else:
        target = "user"
        voice_mute = await database_factory.select(
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            target=target,
            singular=True,
        )
        if voice_mute:
            unvoice_mute_ctx = unvoice_mute_alias_service.UnvoiceMuteMessageContext(
                author_snowflake=message.author.id,
                channel_snowflake=alias_ctx.channel_snowflake,
                guild_snowflake=alias_ctx.guild_snowflake,
                member_snowflake=alias_ctx.member_snowflake,
                message_snowflake=message.id,
                message_channel_snowflake=message.channel.id,
                target=target,
            )
            embed = await unvoice_mute_alias_service.unvoice_mute_by_message(
                ctx=unvoice_mute_ctx, display=True
            )
            return embed
        else:
            voice_mute_ctx = voice_mute_alias_service.VoiceMuteMessageContext(
                author_snowflake=message.author.id,
                channel_snowflake=alias_ctx.channel_snowflake,
                duration_value=alias_ctx.duration_value,
                guild_snowflake=alias_ctx.guild_snowflake,
                member_snowflake=alias_ctx.member_snowflake,
                message_snowflake=message.id,
                message_channel_snowflake=message.channel.id,
                reason=alias_ctx.reason,
                target=target,
            )
            embed = await voice_mute_alias_service.voice_mute_by_message(
                ctx=voice_mute_ctx, display=True
            )
            return embed


async def channel_mute(
    author_snowflake: int,
    channel_snowflake: int,
    duration_value: str,
    guild_snowflake: int,
    reason: str,
    *,
    target: str = "user",
) -> list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    if not isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
        raise commands.CheckFailure("This action must be used in a valid channel.")
    muted_members, pages, skipped_members, failed_members = [], [], [], []
    for member in channel.members:
        if member.id == author_snowflake:
            continue
        try:
            await moderator_service.has_equal_or_lower_role(
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                member_snowflake=author_snowflake,
                target_member_snowflake=member.id,
            )
        except moderator_service.HasEqualOrLowerRole:
            skipped_members.append(member)
            continue
        if member.voice and member.voice.channel:
            if member.voice.channel.id == channel.id:
                try:
                    await member.edit(mute=True)
                except Exception as e:
                    bot.logger.warning(
                        f"Unable to voice-mute member "
                        f"{member.display_name} ({member.id}) in channel "
                        f"{channel.name} ({channel.id}) in guild "
                        f"{channel.guild.name} ({channel.guild.id}). "
                        f"{str(e).capitalize()}"
                    )
                    failed_members.append(member)
        expires_in = duration_builder.parse(duration_value).to_expires_in()
        voice_mute = MODEL(
            channel_snowflake=channel.id,
            expires_in=expires_in,
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
            reason=reason,
            target=target,
        )
        await database_factory.create(voice_mute)
        await voice_mute_alias_service.log_voice_mute(
            author_snowflake=author_snowflake,
            channel_snowflake=channel_snowflake,
            display=True,
            duration_value=duration_value,
            guild_snowflake=guild_snowflake,
            is_channel_scope=True,
            member_snowflake=member.id,
            message_snowflake=None,
            message_channel_snowflake=None,
            reason=reason,
            target=target,
        )
        muted_members.append(member)
    description_lines = [
        f"**Channel:** {channel.mention}",
        f"**Muted:** {len(muted_members)} users",
        f"**Failed:** {len(failed_members)} users",
        f"**Skipped:** {
            len(channel.members) - len(muted_members) - len(failed_members)
        }",
    ]
    embed = discord.Embed(
        description="\n".join(description_lines),
        title=f"{emojis.get_random_emoji()} Channel Mute Summary",
        color=discord.Color.blurple(),
    )
    pages.append(embed)
    return pages


async def channel_unmute(
    channel_snowflake: int,
    guild_snowflake: int,
    *,
    author_snowflake: int | None,
    target: str = "user",
) -> list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    if not isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
        raise commands.CheckFailure("This command must be used in a valid channel.")
    unmuted_members, pages, skipped_members, failed_members = [], [], [], []
    for member in channel.members:
        voice_mute = await database_factory.select(
            channel_snowflake=channel_snowflake,
            member_snowflake=member.id,
            target=target,
            singular=True,
        )
        if not voice_mute:
            skipped_members.append(member)
            continue
        await database_factory.delete(
            target=target,
            channel_snowflake=channel_snowflake,
            member_snowflake=member.id,
        )
        if member.voice and member.voice.channel:
            if member.voice.channel.id == channel.id:
                try:
                    await member.edit(mute=False)
                except Exception as e:
                    bot.logger.warning(
                        f"Unable to undo voice-mute "
                        f"for member {member.display_name} ({member.id}) "
                        f"in channel {channel.name} ({channel.id}) "
                        f"in guild {channel.guild.name} "
                        f"({channel.guild.id}). "
                        f"{str(e).capitalize()}"
                    )
                    failed_members.append(member)
        await unvoice_mute_alias_service.log_unvoice_mute(
            author_snowflake=author_snowflake,
            channel_snowflake=channel_snowflake,
            display=True,
            guild_snowflake=guild_snowflake,
            is_channel_scope=True,
            member_snowflake=member.id,
            message_snowflake=None,
            message_channel_snowflake=None,
            target=target,
        )
        unmuted_members.append(member)
    description_lines = [
        f"**Channel:** {channel.mention}",
        f"**Unmuted:** {len(unmuted_members)} users",
        f"**Failed:** {len(failed_members)} users",
        f"**Skipped:** {
            len(channel.members) - len(unmuted_members) - len(failed_members)
        }",
    ]
    embed = discord.Embed(
        description="\n".join(description_lines),
        title=f"{emojis.get_random_emoji()} Channel Unmute Summary",
        color=discord.Color.blurple(),
    )
    pages.append(embed)
    return pages


async def is_voice_muted(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    voice_mute = await database_factory.select(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        singular=True,
    )
    if voice_mute:
        return True
    return False
