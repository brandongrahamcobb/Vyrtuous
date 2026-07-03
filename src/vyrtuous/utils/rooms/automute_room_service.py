"""!/bin/python3stage"
automute_room_service.py The purpose of this program is to extend Service to service the stage class.

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

import time
from dataclasses import dataclass, field
from typing import Any, Dict, List

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.automute import Stage
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.moderation import voice_mute_service
from vyrtuous.utils.users import moderator_service

MODEL = Stage


@dataclass
class StageDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, Any]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


async def send_stage_ask_to_speak_message(
    join_log: dict[int, discord.Member], member: discord.Member, stage: Stage
):
    bot = DiscordBot.get_instance()
    now = time.time()
    join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
    if len(join_log[member.id]) < 1:
        join_log[member.id].append(now)
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} {stage.channel_snowflake} — Stage Mode",
            description=f"Ends <t:{int(stage.expires_in.timestamp())}:R>",
            color=discord.Color.green(),
        )
        embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
        await bot.get_channel(stage.channel_snowflake).send(embed=embed)


async def build_dictionary(obj):
    database_factory = DatabaseFactory(MODEL)
    stages = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        stages = await database_factory.select(guild_snowflake=obj.id, singular=False)
    elif isinstance(obj, discord.abc.GuildChannel):
        stages = await database_factory.select(channel_snowflake=obj.id, singular=False)
    else:
        stages = await database_factory.select(singular=False)
    if stages:
        for stage in stages:
            dictionary.setdefault(stage.guild_snowflake, {"channels": {}})
            dictionary[stage.guild_snowflake]["channels"].setdefault(
                stage.channel_snowflake, {}
            )
            dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ].setdefault("stages", {})
            dictionary[stage.guild_snowflake]["channels"][stage.channel_snowflake][
                "stages"
            ].update({"expires_in": stage.expires_in})
    return dictionary


async def build_pages(is_at_home: bool, obj):
    bot = DiscordBot.get_instance()
    lines, pages = [], []

    obj_name = "All Servers"
    if obj:
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Stages for {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary = await list_service.process_dictionary(
        cls=StageDictionary, dictionary=dictionary
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        stage_n = 0
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for channel_snowflake, stage_dictionary in guild_data.get("channels").items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            lines.append(
                f"**Expires in:** {stage_dictionary.get('stages', {}).get('expires_in', None)}"
            )
            stage_n += 1
            field_count += 1
            if field_count == list_service.CHUNK_SIZE:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed = list_service.flush_page(embed, pages, title, guild.name)
                lines = []
                field_count = 0
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({stage_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_channels)
        pages.extend(processed_dictionary.skipped_guilds)
    if not pages:
        return "No stages found."
    return pages


async def toggle_stage(channel, context, duration_value):
    database_factory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    failed, pages, skipped, succeeded = [], [], [], []
    stage = await database_factory.select(channel_snowflake=channel.id, singular=True)
    if stage:
        title = f"{emojis.get_random_emoji()} Stage Ended in {channel.mention}"
        await database_factory.delete(channel_snowflake=channel.id)
        failed, succeeded = await voice_mute_service.off_stage(channel=channel)
        description_lines = [
            f"**Channel:** {channel.mention}",
            f"**Unmuted:** {len(succeeded)} users",
        ]
        if failed:
            description_lines.append(f"**Failed:** {len(failed)}")
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=title,
            color=discord.Color.blurple(),
        )
        pages.append(embed)
    else:
        stage = MODEL(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            expires_in=duration_builder.parse(value=duration_value).to_expires_in(),
        )
        await database_factory.create(stage)
        failed, skipped, succeeded = await voice_mute_service.on_stage(
            channel=channel,
            context=context,
            duration_value=duration_value,
        )
        description_lines = [
            f"**Channel:** {channel.mention}",
            f"**Expires:** {duration_builder.parse(value=duration_value).to_unix_ts()}",
            f"**Muted:** {len(succeeded)} users",
            f"**Skipped:** {len(skipped)}",
        ]
        if failed:
            description_lines.append(f"**Failed:** {len(failed)}")
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{emojis.get_random_emoji()} Stage Created in {channel.name}",
            color=discord.Color.blurple(),
        )
        pages.append(embed)
    return pages


async def toggle_stage_mute(channel, context, member):
    database_factory = DatabaseFactory(MODEL)
    await moderator_service.has_equal_or_lower_role(
        **context.to_dict(),
        target_member_snowflake=member.id,
    )
    stage = await database_factory.select(
        channel_snowflake=channel.id,
        guild_snowflake=channel.guild.id,
        singular=True,
    )
    if stage:
        await member.edit(mute=not member.voice.mute)
        return (
            f"Successfully toggled the mute for {member.mention} in {channel.mention}."
        )


async def clean_expired():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    expired_stages = await database_factory.select(expired=True, singular=False)
    if expired_stages:
        for expired_stage in expired_stages:
            channel_snowflake = int(expired_stage.channel_snowflake)
            guild_snowflake = int(expired_stage.guild_snowflake)
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                await database_factory.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                )
                bot.logger.info(
                    f"Unable to locate guild {guild_snowflake}, cleaning up expired stage."
                )
                continue
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                await database_factory.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                )
                bot.logger.info(
                    f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                )
                continue
            voice_mute_service.clean_stage_expired(channel=channel, guild=guild)
        bot.logger.info("Cleaned up expired stages.")


async def migrate(kwargs):
    database_factory = DatabaseFactory(MODEL)
    await database_factory.update(**kwargs)


async def enforce(after, member):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    should_be_muted = False
    expires_in = None
    stage = await database_factory.select(
        channel_snowlfake=after.channel, singular=True
    )
    if stage:
        await send_stage_ask_to_speak_message(
            join_log=bot.join_log, member=member, stage=stage
        )
        highest_role = await moderator_service.resolve_highest_role(
            channel_snowflake=after.channel.id,
            guild_snowflake=after.channel.guild.id,
            member_snowflake=member.id,
        )
        if highest_role == "Everyone":
            should_be_muted = True
            expires_in = stage.expires_in
    return should_be_muted, expires_in


async def is_active_stage_room(channel):
    database_factory = DatabaseFactory(MODEL)
    stage_room = await database_factory.select(
        channel_snowflake=channel.id, singular=True
    )
    if stage_room:
        return True
    return False
