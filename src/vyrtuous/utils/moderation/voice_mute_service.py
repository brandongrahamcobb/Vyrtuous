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

from datetime import datetime, timedelta, timezone

import discord

from vyrtuous.aliases import unvoice_mute_alias_service, voice_mute_alias_service
from vyrtuous.aliases.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot
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
    target = "user"
    database_factory = DatabaseFactory(MODEL)
    voice_mute = await database_factory.select(
        channel_snowflake=alias_ctx.channel_snowflake,
        guild_snowflake=alias_ctx.guild_snowflake,
        member_snowflake=alias_ctx.member_snowflake,
        target=target,
        singular=True,
    )
    if voice_mute:
        ctx = unvoice_mute_alias_service.UnvoiceMuteMessageContext(
            author_snowflake=message.author.id,
            channel_snowflake=alias_ctx.channel_snowflake,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            target=target,
        )
        embed = await unvoice_mute_alias_service.unvoice_mute_by_message(
            ctx=ctx, display=True
        )
        return embed
    else:
        ctx = voice_mute_alias_service.VoiceMuteMessageContext(
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
            ctx=ctx, display=True
        )
        return embed


async def clean_expired():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    expired_voice_mutes = await database_factory.select(expired=True, singular=False)
    if expired_voice_mutes:
        for expired_voice_mute in expired_voice_mutes:
            channel_snowflake = int(expired_voice_mute.channel_snowflake)
            guild_snowflake = int(expired_voice_mute.guild_snowflake)
            member_snowflake = int(expired_voice_mute.member_snowflake)
            target = expired_voice_mute.target
            guild = bot.get_guild(guild_snowflake)
            kwargs = {
                "channel_snowflake": channel_snowflake,
                "guild_snowflake": guild_snowflake,
                "member_snowflake": member_snowflake,
                "target": target,
            }
            if guild is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate guild {guild_snowflake}, cleaning up expired voice-mute."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                )
                continue
            member = guild.get_member(member_snowflake)
            if member is None:
                await database_factory.delete(**kwargs)
                bot.logger.info(
                    f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.name}), cleaning up expired voice-mute."
                )
                continue
            await database_factory.delete(**kwargs)
            if (
                member.voice
                and member.voice.channel
                and member.voice.channel.id == channel_snowflake
            ):
                try:
                    await member.edit(mute=False)
                    bot.logger.info(
                        f"Undone voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake})."
                    )
                except discord.Forbidden as e:
                    bot.logger.warning(
                        f"Unable to undo voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
                    )
            else:
                bot.logger.info(
                    f"Member {member.display_name} ({member.id}) is not in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), skipping undo voice-mute."
                )


async def room_mute(author, channel, reason):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    muted_members, pages, skipped_members, failed_members = [], [], [], []
    for member in channel.members:
        if member.id == author.id:
            continue
        voice_mute = await database_factory.select(
            channel_snowflake=channel.id,
            member_snowflake=member.id,
            target="user",
            singular=True,
        )
        try:
            await moderator_service.has_equal_or_lower_role(
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                member_snowflake=author.id,
                target_member_snowflake=member.id,
            )
        except moderator_service.HasEqualOrLowerRole:
            skipped_members.append(member)
            continue
        if voice_mute:
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
        expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
        voice_mute = MODEL(
            channel_snowflake=channel.id,
            expires_in=expires_in,
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
            reason=reason,
            target="user",
        )
        await database_factory.create(voice_mute)
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
        title=f"{emojis.get_random_emoji()} Room Mute Summary",
        color=discord.Color.blurple(),
    )
    pages.append(embed)
    return pages


async def room_unmute(channel):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    unmuted_members, pages, skipped_members, failed_members = [], [], [], []
    for member in channel.members:
        voice_mute = await database_factory.select(
            channel_snowflake=channel.id,
            member_snowflake=member.id,
            target="user",
            singular=True,
        )
        if not voice_mute:
            skipped_members.append(member)
            continue
        await database_factory.delete(
            target="user", channel_snowflake=channel.id, member_snowflake=member.id
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
        title=f"{emojis.get_random_emoji()} Room Unmute Summary",
        color=discord.Color.blurple(),
    )
    pages.append(embed)
    return pages


async def clean_expired_stage(channel, guild):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    voice_mutes = await database_factory.select(
        channel_snowflake=channel.id,
        guild_snowflake=guild.id,
        target="room",
        singular=False,
    )
    for voice_mute in voice_mutes:
        member_snowflake = voice_mute.member_snowflake
        member = guild.get_member(member_snowflake)
        if member is None:
            await database_factory.delete(
                channel_snowflake=channel.id,
                member_snowflake=int(member_snowflake),
                guild_snowflake=guild.id,
                target="room",
            )
            bot.logger.info(
                f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.id}) from expired stage."
            )
            continue
        await database_factory.delete(
            channel_snowflake=channel.id,
            member_snowflake=member.id,
            guild_snowflake=guild.id,
            target="room",
        )
        if (
            member.voice
            and member.voice.channel
            and member.voice.mute
            and member.voice_channel.id == channel.id
        ):
            try:
                await member.edit(mute=False, reason="Stage room closed automatically.")
                bot.logger.info(
                    f"Undone voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in in guild {guild.name} ({guild.id}) after stage expired."
                )
            except discord.Forbidden as e:
                bot.logger.warning(
                    f"Unable to undo voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.id}). {str(e).capitalize()}"
                )
        else:
            bot.logger.info(
                f"Member {member.display_name} ({member.id}) is not in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.id}), skipping undo voice-mute."
            )


async def off_stage(channel):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    failed, succeeded = [], []
    for member in channel.members:
        await database_factory.delete(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
            target="room",
        )
        voice_mute = await database_factory.select(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
            target="user",
            singular=True,
        )
        if not voice_mute and member.voice and member.voice.mute:
            try:
                await member.edit(
                    mute=False,
                    reason="Stage ended — no user-specific mute found",
                )
                succeeded.append(member)
            except discord.Forbidden as e:
                bot.logger.warning(
                    f"Unable to undo voice-mute "
                    f"for member {member.display_name} ({member.id}) in "
                    f"channel {channel.name} ({channel.id}) in "
                    f"guild {channel.guild.name} ({channel.guild.id}). "
                    f"{str(e).capitalize()}"
                )
                failed.append(member)
    return failed, succeeded


async def on_stage(channel, context, duration_value):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    failed, skipped, succeeded = [], [], []
    for member in channel.members:
        if (
            await moderator_service.check_minimum_role(
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                member_snowflake=member.id,
                lowest_role="Coordinator",
            )
            or member.id == context.author.id
        ):
            skipped.append(member)
            continue
        voice_mute = MODEL(
            channel_snowflake=channel.id,
            expires_in=duration_builder.parse(value=duration_value).to_expires_in(),
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
            target="room",
            reason="Stage mute",
        )
        await database_factory.create(voice_mute)
        try:
            if member.voice and member.voice.channel.id == channel.id:
                await member.edit(mute=True)
            succeeded.append(member)
        except Exception as e:
            bot.logger.warning(
                f"Unable to voice-mute "
                f"member {member.display_name} ({member.id}) "
                f"in channel {channel.name} ({channel.id}) "
                f"in guild {channel.guild.name} ({channel.guild.id}). "
                f"{str(e).capitalize()}"
            )
            failed.append(member)
    return failed, skipped, succeeded


async def is_voice_muted(channel, member):
    database_factory = DatabaseFactory(MODEL)
    voice_mute = await database_factory.select(
        channel_snowflake=channel.id, member_snowflake=member.id, singular=True
    )
    if voice_mute:
        return True
    return False
