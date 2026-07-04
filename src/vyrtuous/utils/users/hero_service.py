"""!/bin/python3
permission_service.py The purpose of this program is to provide the service for deciding whether a member has sufficient permissions.

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
from vyrtuous.cache.hero import Hero
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)

INFRACTION_MODELS = [
    ban_service.MODEL,
    flag_service.MODEL,
    text_mute_service.MODEL,
    voice_mute_service.MODEL,
]

MODEL = Hero


async def unrestrict(guild_snowflake: int, member_snowflake: int) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if not guild:
        raise commands.GuildNotFound(str(guild_snowflake))
    member = guild.get_member(member_snowflake)
    if not member:
        raise commands.MemberNotFound(str(member_snowflake))
    kwargs = {
        "guild_snowflake": int(guild_snowflake),
        "member_snowflake": member_snowflake,
    }
    for model in INFRACTION_MODELS:
        database_factory: DatabaseFactory = DatabaseFactory(model)
        objects = await database_factory.select(singular=False, **kwargs)
        if objects:
            match model.identifier:
                case "ban":
                    for ban in objects:
                        channel = guild.get_channel(ban.channel_snowflake)
                        if channel:
                            try:
                                await channel.set_permissions(member, overwrite=None)
                            except discord.Forbidden:
                                bot.logger.warning(
                                    f"Unable to unban {member.name} ({member.id}) in {channel.name} ({channel.id})."
                                )
                case "tmute":
                    for text_mute in objects:
                        channel = guild.get_channel(text_mute.channel_snowflake)
                        if channel:
                            try:
                                await channel.set_permissions(
                                    member, send_messages=True
                                )
                            except discord.Forbidden:
                                bot.logger.warning(
                                    f"Unable to untmute {member.name} ({member.id}) in {channel.name} ({channel.id})."
                                )
                case "vmute":
                    for voice_mute in objects:
                        channel = guild.get_channel(voice_mute.channel_snowflake)
                        if channel and member.voice and member.voice.mute:
                            try:
                                await member.edit(mute=False)
                            except discord.Forbidden:
                                bot.logger.warning(
                                    f"Unable to unmute {member.name} ({member.id}) in {channel.name} ({channel.id})."
                                )
            await database_factory.delete(**kwargs)


async def add_invincible_member(guild_snowflake: int, member_snowflake: int) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    hero = Hero(
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
    await database_factory.create(hero)
    bot.registry.get(MemberState).invincible[guild_snowflake].add(member_snowflake)


async def remove_invincible_member(guild_snowflake: int, member_snowflake: int) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    await database_factory.delete(
        guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
    )
    bot.registry.get(MemberState).invincible[guild_snowflake].discard(member_snowflake)
