"""!/bin/python3
voice_mute_service.py The purpose of this program is to service voice-mutes.

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

from vyrtuous.aliases import unvoice_mute_alias_service, voice_mute_alias_service
from vyrtuous.aliases.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState, MemberState, PermissionState
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.utils.errors.error import (
    ChannelNotFound,
    CheckFailure,
    GuildNotFound,
    HasEqualOrLowerRole,
)
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.permissions import permission_service

MODEL = VoiceMute


async def enforce_or_undo(
    alias_ctx: AliasContext,
    message: discord.Message,
) -> discord.Embed:
    duration_builder = DurationBuilder()
    automute_database_factory: DatabaseFactory = DatabaseFactory(AutoMute)
    auto_mute_channel = await automute_database_factory.select(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        singular=True,
    )
    if await is_voice_muted(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        member_snowflake=alias_ctx.member_snowflake,
        targets=["command"],
    ):
        unvoice_mute_ctx = unvoice_mute_alias_service.UnvoiceMuteMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            target="command",
        )
        embed = await unvoice_mute_alias_service.unvoice_mute_by_message(
            ctx=unvoice_mute_ctx, display=True
        )
        return embed
    elif auto_mute_channel:
        target = "auto"
        voice_mute_ctx = voice_mute_alias_service.VoiceMuteMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            duration=duration_builder.from_timestamp(
                auto_mute_channel.expires_in
            ).build(),
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            reason="Automuted",
            target=target,
        )
        embed = await voice_mute_alias_service.voice_mute_by_message(
            ctx=voice_mute_ctx, display=True
        )
        return embed
    else:
        target = "command"
        voice_mute_ctx = voice_mute_alias_service.VoiceMuteMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            duration=alias_ctx.duration,
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
    duration: DurationObject,
    excluded: list[int],
    guild_snowflake: int,
    reason: str,
    *,
    target: str = "click",
) -> list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise ChannelNotFound(str(channel_snowflake))
    if not isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
        raise CheckFailure("This command must target a valid channel.")
    muted_members, pages, skipped_members, failed_members = [], [], [], []
    for member in channel.members:
        if member.id in excluded:
            continue
        try:
            await permission_service.has_equal_or_lower_role(
                permission_state=permission_state,
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                author_snowflake=author_snowflake,
                member_snowflake=member.id,
            )
        except HasEqualOrLowerRole:
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
        expires_in = duration_builder.load(duration=duration).to_expires_in()
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
            duration=duration,
            guild_snowflake=guild_snowflake,
            is_channel_scope=True,
            member_snowflake=member.id,
            message_snowflake=None,
            message_channel_snowflake=None,
            reason=reason,
            target=target,
        )
        if target == "auto":
            bot.registry.get(MemberState).automuted.setdefault(member.id, set()).add(
                channel_snowflake
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
    author_snowflake: int,
    channel_snowflake: int,
    guild_snowflake: int,
    reason: str | None = "No reason provided.",
    target: str = "click",
) -> list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    permission_state: PermissionState = bot.registry.get(PermissionState)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise ChannelNotFound(str(channel_snowflake))
    if not isinstance(channel, (discord.VoiceChannel, discord.StageChannel)):
        raise CheckFailure("This command must target a valid channel.")
    unmuted_members, pages, skipped_members, failed_members = [], [], [], []
    for member in channel.members:
        try:
            await permission_service.has_equal_or_lower_role(
                permission_state=permission_state,
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                author_snowflake=author_snowflake,
                member_snowflake=member.id,
            )
        except HasEqualOrLowerRole:
            skipped_members.append(member)
            continue
        voice_mutes = await database_factory.select(
            channel_snowflake=channel_snowflake,
            member_snowflake=member.id,
            singular=True,
        )
        await database_factory.delete(
            target=target,
            channel_snowflake=channel_snowflake,
            member_snowflake=member.id,
        )
        if not voice_mutes:
            skipped_members.append(member)
            continue
        if len(voice_mutes) and member.voice and member.voice.channel:
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
            reason=reason or "No reason provided.",
            target=target,
        )
        if target == "auto":
            member_dict = bot.registry.get(MemberState).automuted
            for channel_snowflakes in member_dict.values():
                channel_snowflakes.discard(channel_snowflake)
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
    guild_snowflake: int,
    member_snowflake: int,
    targets: list[str],
    channel_snowflake: int | None = None,
) -> list[str]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    muted_targets = []
    for target in targets:
        if channel_snowflake is None:
            voice_mute = await database_factory.select(
                guild_snowflake=guild_snowflake,
                member_snowflake=member_snowflake,
                target=target,
                singular=True,
            )
        else:
            voice_mute = await database_factory.select(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                member_snowflake=member_snowflake,
                target=target,
                singular=True,
            )
        if voice_mute:
            muted_targets.append(target)
    return muted_targets


async def populate(target: str = "auto"):
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    voice_mutes = await database_factory.select(target=target, singular=False)
    for voice_mute in voice_mutes:
        if target == "auto":
            bot.registry.get(MemberState).automuted.setdefault(
                voice_mute.member_snowflake, set()
            ).add(voice_mute.channel_snowflake)


async def alert_mute(
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
    target: str = "command",
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
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    voice_mute = await database_factory.select(
        channel_snowflake=channel_snowflake,
        member_snowflake=member_snowflake,
        target=target,
        singular=True,
    )
    if voice_mute and bot.registry.get(ChannelState).should_notify(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        timeout=900.0,
    ):
        if isinstance(channel, discord.channel.VocalGuildChannel):
            duration_builder = DurationBuilder()
            embed = discord.Embed(
                title=f"\u274c {member.display_name} must be unmuted via /mute.",
                description=f"**Expires:** {duration_builder.from_timestamp(voice_mute.expires_in).to_unix_ts()}\n**Reason:** {voice_mute.reason}",
                color=discord.Color.red(),
            )
            embed.set_thumbnail(url=member.display_avatar.url)
            await channel.send(embed=embed)
    bot.registry.get(ChannelState).record(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
