"""!/bin/python3stage"
automute_channel_service.py The purpose of this program is to extend Service to service the stage class.

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


# async def send_automute_ask_to_speak_message(
#     join_log: dict[int, discord.Member], member: discord.Member, automute: AutoMute
# ):
#     bot: DiscordBot = DiscordBot.get_instance()
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

# async def toggle_automute(channel_snowflake: int, guild_snowflake: int, context, duration_value):
#     bot: DiscordBot = DiscordBot.get_instance()
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


# async def enforce(after, member):
#     bot: DiscordBot = DiscordBot.get_instance()
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


async def is_active_automute_channel(channel):
    database_factory = DatabaseFactory(MODEL)
    automute_channel = await database_factory.select(
        channel_snowflake=channel.id, singular=True
    )
    if automute_channel:
        return True
    return False
