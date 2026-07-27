"""!/bin/python3
cap_service.py The purpose of this program is to extend Service to service the cap command class.

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

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.cap import Cap
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.models.duration import DurationBuilder, DurationObject

MODEL = Cap


async def toggle_cap(
    category: str,
    channel_snowflake: int,
    duration: DurationObject,
    guild_snowflake: int,
) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(str(channel_snowflake))
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    duration_builder: DurationBuilder = DurationBuilder()
    seconds = duration_builder.load(duration).to_seconds()
    where_kwargs = {
        "channel_snowflake": channel_snowflake,
        "category": category,
        "guild_snowflake": guild_snowflake,
    }
    cap = await database_factory.select(
        singular=True, channel_snowflake=channel_snowflake, category=category
    )
    if cap and seconds:
        await database_factory.update(
            set_kwargs={"duration_seconds": seconds}, where_kwargs=where_kwargs
        )
        return f"Cap `{category}` modified for {channel.mention}."
    elif cap:
        await database_factory.delete(
            channel_snowflake=channel_snowflake, category=category
        )
        return (
            f"Cap of type {category} "
            f"and channel {channel.mention} deleted successfully."
        )
    else:
        where_kwargs.update({"duration_seconds": seconds})
        cap = MODEL(
            channel_snowflake=channel.id,
            category=category,
            duration_seconds=seconds,
            guild_snowflake=guild_snowflake,
        )
        await database_factory.create(cap)
        return f"Cap `{category}` created for {channel.mention} successfully."


async def exceeds_cap(
    category: str,
    channel_snowflake: int,
    duration: DurationObject,
    guild_snowflake: int,
) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    duration_builder = DurationBuilder()
    duration_seconds = duration_builder.load(duration=duration).to_seconds()
    exceeds_cap = False
    cap = await database_factory.select(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        category=category,
        singular=True,
    )
    number = duration.number
    if cap:
        if duration_seconds > cap.duration_seconds or number == 0:
            exceeds_cap = True
    else:
        cap_duration_seconds = duration_builder.parse(value="8h").to_seconds()
        if duration_seconds > cap_duration_seconds or number == 0:
            exceeds_cap = True
    return exceeds_cap
