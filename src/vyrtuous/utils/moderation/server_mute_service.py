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
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.utils.users import moderator_service

MODEL = VoiceMute


async def toggle_server_mute(
    author_snowflake: int,
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
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member:
            member_name = simplified_member[0]
        else:
            raise commands.MemberNotFound(str(member_snowflake))
    else:
        member_name = member.mention
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    await moderator_service.has_equal_or_lower_role(
        channel_snowflake=None,
        guild_snowflake=guild_snowflake,
        member_snowflake=author_snowflake,
        target_member_snowflake=member_snowflake,
    )
    server_mute = await database_factory.select(
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        target="server",
        singular=True,
    )
    if not server_mute:
        server_mute = VoiceMute(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            reason=reason,
            target="server",
        )
        await database_factory.create(server_mute)
        action = "muted"
        should_be_muted = True
    else:
        await database_factory.delete(
            guild_snowflake=guild_snowflake,
            member_snowflake=member_snowflake,
            target="server",
        )
        action = "unmuted"
        should_be_muted = False
    if member and member.voice and member.voice.channel:
        await member.edit(mute=should_be_muted)
    return f"Successfully server {action} {member_name} in {guild.name}."
