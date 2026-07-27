"""!/bin/python3stage"
automute_channel_service.py The purpose of this program is to extend Service to service the stage class.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import Duration, DurationBuilder, DurationObject
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.moderation import voice_mute_service

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


async def toggle_automute(
    author_snowflake: int,
    channel_snowflake: int,
    duration: DurationObject | None,
    guild_snowflake: int,
):
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    pages: list[discord.Embed] = []
    automute = await database_factory.select(
        channel_snowflake=channel_snowflake, singular=True
    )
    if automute:
        await database_factory.delete(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
        )
        embed = await voice_mute_service.channel_unmute(
            author_snowflake=author_snowflake,
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            target="auto",
        )
        pages.extend(embed)
    else:
        if duration:
            automute = MODEL(
                channel_snowflake=channel_snowflake,
                guild_snowflake=guild_snowflake,
                expires_in=duration_builder.load(duration=duration).to_expires_in(),
            )
            await database_factory.create(automute)
            embed = await voice_mute_service.channel_mute(
                author_snowflake=author_snowflake,
                channel_snowflake=channel_snowflake,
                duration=duration,
                excluded=[author_snowflake],
                guild_snowflake=guild_snowflake,
                reason=f"Automute enabled.",
                target="auto",
            )
            pages.extend(embed)
    return pages


async def is_active_automute_channel(channel_snowflake: int, guild_snowflake: int):
    database_factory = DatabaseFactory(MODEL)
    automute_channel = await database_factory.select(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        singular=True,
    )
    if automute_channel:
        return True
    return False
