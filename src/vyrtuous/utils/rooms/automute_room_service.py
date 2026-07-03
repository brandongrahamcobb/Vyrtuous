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

from dataclasses import dataclass, field
from typing import Any, Dict, List

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.listing import list_service
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.moderation import voice_mute_service
from vyrtuous.utils.users import moderator_service

MODEL = AutoMute


@dataclass
class AutoMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, Any]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


# async def send_automute_ask_to_speak_message(
#     join_log: dict[int, discord.Member], member: discord.Member, automute: AutoMute
# ):
#     bot = DiscordBot.get_instance()
#     now = time.time()
#     join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
#     if len(join_log[member.id]) < 1:
#         join_log[member.id].append(now)
#         embed = discord.Embed(
#             title=f"{emojis.get_random_emoji()} {automute.channel_snowflake} — Stage Mode",
#             description=f"Ends <t:{int(automute.expires_in.timestamp())}:R>",
#             color=discord.Color.green(),
#         )
#         embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
#         await bot.get_channel(automute.channel_snowflake).send(embed=embed)


async def build_dictionary(obj):
    database_factory = DatabaseFactory(MODEL)
    automutes = []
    dictionary = {}
    if isinstance(obj, discord.Guild):
        automutes = await database_factory.select(
            guild_snowflake=obj.id, singular=False
        )
    elif isinstance(obj, discord.abc.GuildChannel):
        automutes = await database_factory.select(
            channel_snowflake=obj.id, singular=False
        )
    else:
        automutes = await database_factory.select(singular=False)
    if automutes:
        for automute in automutes:
            dictionary.setdefault(automute.guild_snowflake, {"channels": {}})
            dictionary[automute.guild_snowflake]["channels"].setdefault(
                automute.channel_snowflake, {}
            )
            dictionary[automute.guild_snowflake]["channels"][
                automute.channel_snowflake
            ].setdefault("automutes", {})
            dictionary[automute.guild_snowflake]["channels"][
                automute.channel_snowflake
            ]["automutes"].update({"expires_in": automute.expires_in})
    return dictionary


async def build_pages(is_at_home: bool, obj):
    bot = DiscordBot.get_instance()
    lines, pages = [], []

    obj_name = "All Servers"
    if obj:
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Automute Rooms for {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary = await list_service.process_dictionary(
        cls=AutoMuteDictionary, dictionary=dictionary
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        automute_n = 0
        field_count = 0
        lines = []
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for channel_snowflake, automute_dictionary in guild_data.get(
            "channels"
        ).items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            lines.append(
                f"**Expires in:** {automute_dictionary.get('automutes', {}).get('expires_in', None)}"
            )
            automute_n += 1
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
        embed.description = f"**{original_description} ({automute_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_channels)
        pages.extend(processed_dictionary.skipped_guilds)
    if not pages:
        return "No automute rooms found."
    return pages


# async def toggle_automute(channel_snowflake: int, guild_snowflake: int, context, duration_value):
#     bot = DiscordBot.get_instance()
#     guild = bot.get_guild(guild_snowflake)
#     if guild is None:
#         raise commands.GuildNotFound(str(guild_snowflake))
#     channel = guild.get_channel(channel_snowflake)
#     if channel is None:
#         raise commands.ChannelNotFound(str(channel_snowflake))
#     database_factory = DatabaseFactory(MODEL)
#     duration_builder = DurationBuilder()
#     failed, pages, skipped, succeeded = [], [], [], []
#     automute = await database_factory.select(
#         channel_snowflake=channel_snowflake, singular=True
#     )
#     if automute:
#         title = f"{emojis.get_random_emoji()} Stage Ended in {channel.mention}"
#         await database_factory.delete(channel_snowflake=channel.id)
#         for
#         failed, succeeded = await voice_mute_service.off_automute(channel=channel)
#         description_lines = [
#             f"**Channel:** {channel.mention}",
#             f"**Unmuted:** {len(succeeded)} users",
#         ]
#         if failed:
#             description_lines.append(f"**Failed:** {len(failed)}")
#         embed = discord.Embed(
#             description="\n".join(description_lines),
#             title=title,
#             color=discord.Color.blurple(),
#         )
#         pages.append(embed)
#     else:
#         automute = MODEL(
#             channel_snowflake=channel.id,
#             guild_snowflake=channel.guild.id,
#             expires_in=duration_builder.parse(value=duration_value).to_expires_in(),
#         )
#         await database_factory.create(automute)
#         failed, skipped, succeeded = await voice_mute_service.on_stage(
#             channel=channel,
#             context=context,
#             duration_value=duration_value,
#         )
#         description_lines = [
#             f"**Channel:** {channel.mention}",
#             f"**Expires:** {duration_builder.parse(value=duration_value).to_unix_ts()}",
#             f"**Muted:** {len(succeeded)} users",
#             f"**Skipped:** {len(skipped)}",
#         ]
#         if failed:
#             description_lines.append(f"**Failed:** {len(failed)}")
#         embed = discord.Embed(
#             description="\n".join(description_lines),
#             title=f"{emojis.get_random_emoji()} Stage Created in {channel.name}",
#             color=discord.Color.blurple(),
#         )
#         pages.append(embed)
#     return pages
#
#
# async def toggle_automute_mute(channel, context, member):
#     database_factory = DatabaseFactory(MODEL)
#     await moderator_service.has_equal_or_lower_role(
#         **context.to_dict(),
#         target_member_snowflake=member.id,
#     )
#     automute = await database_factory.select(
#         channel_snowflake=channel.id,
#         guild_snowflake=channel.guild.id,
#         singular=True,
#     )
#     if automute:
#         await member.edit(mute=not member.voice.mute)
#         return (
#             f"Successfully toggled the mute for {member.mention} in {channel.mention}."
#         )


async def clean_expired():
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    expired_automutes = await database_factory.select(expired=True, singular=False)
    if expired_automutes:
        for expired_automute in expired_automutes:
            channel_snowflake = int(expired_automute.channel_snowflake)
            guild_snowflake = int(expired_automute.guild_snowflake)
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                await database_factory.delete(
                    channel_snowflake=channel_snowflake,
                    guild_snowflake=guild_snowflake,
                )
                bot.logger.info(
                    f"Unable to locate guild {guild_snowflake}, cleaning up expired automute."
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
            if not isinstance(channel, discord.VoiceChannel):
                continue
            database_factory = DatabaseFactory(VoiceMute)
            automutes = await database_factory.select(
                channel_snowflake=channel_snowflake, target="auto", singular=False
            )
            for automute in automutes:
                member_snowflake = automute.member_snowflake
                member = guild.get_member(member_snowflake)
                if member is None:
                    continue
                else:
                    await database_factory.delete(
                        channel_snowflake=channel_snowflake,
                        member_snowflake=member_snowflake,
                        target="auto",
                    )
                    if member in channel.members:
                        try:
                            await member.edit(mute=False, reason="Undoing automute")
                        except discord.Forbidden:
                            continue
        bot.logger.info("Cleaned up expired automutes.")


# async def enforce(after, member):
#     bot = DiscordBot.get_instance()
#     database_factory = DatabaseFactory(MODEL)
#     should_be_muted = False
#     expires_in = None
#     automute = await database_factory.select(
#         channel_snowlfake=after.channel, singular=True
#     )
#     if automute:
#         await send_automute_ask_to_speak_message(
#             join_log=bot.join_log, member=member, automute=automute
#         )
#         highest_role = await moderator_service.resolve_highest_role(
#             channel_snowflake=after.channel.id,
#             guild_snowflake=after.channel.guild.id,
#             member_snowflake=member.id,
#         )
#         if highest_role == "Everyone":
#             should_be_muted = True
#             expires_in = automute.expires_in
#     return should_be_muted, expires_in


async def is_active_automute_room(channel):
    database_factory = DatabaseFactory(MODEL)
    automute_room = await database_factory.select(
        channel_snowflake=channel.id, singular=True
    )
    if automute_room:
        return True
    return False
