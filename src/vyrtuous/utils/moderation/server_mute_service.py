"""!/bin/python3
server_mute_service.py The purpose of this program is to extend AliasService to service the server mute infraction.

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

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.utils.users import moderator_service

MODEL = VoiceMute


async def toggle_server_mute(
    author_snowflake: int,
    channel_snowflake: int,
    guild_snowflake: int,
    member_snowflake: int,
    reason: str,
) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(guild_snowflake))
    member = guild.get_member(member_snowflake)
    if member is None:
        raise commands.MemberNotFound(str(member_snowflake))
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    await moderator_service.has_equal_or_lower_role(
        channel_snowflake=channel_snowflake,
        guild_snowflake=guild_snowflake,
        member_snowflake=author_snowflake,
        target_member_snowflake=member_snowflake,
    )
    server_mute = await database_factory.select(
        singular=True,
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
    )
    if not server_mute:
        server_mute = VoiceMute(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            reason=reason,
        )
        await database_factory.create(server_mute)
        action = "muted"
        should_be_muted = True
    else:
        await database_factory.delete(
            guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
        )
        action = "unmuted"
        should_be_muted = False
    if member.voice and member.voice.channel:
        await member.edit(mute=should_be_muted)
    return f"Successfully server {action} {member.mention} in {guild.name}."


async def is_server_muted(channel, member) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    server_mute = await database_factory.select(
        channel_snowflake=channel.id, member_snowflake=member.id, singular=False
    )
    if server_mute:
        return True
    return False
