"""!/bin/python3stage"
automute_channel_service.py The purpose of this program is service automute channels.

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

from vyrtuous.db.automute import AutoMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.utils.moderation import voice_mute_service
from vyrtuous.utils.tracking import data_builder, stream_service

MODEL = AutoMute


async def toggle_automute(
    author_snowflake: int,
    channel_snowflake: int,
    duration: DurationObject | None,
    guild_snowflake: int,
):
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
