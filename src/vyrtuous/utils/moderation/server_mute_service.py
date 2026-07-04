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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.server_mute import ServerMute
from vyrtuous.utils.users import moderator_service

MODEL = ServerMute


async def toggle_server_mute(context, member, reason) -> str:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    await moderator_service.has_equal_or_lower_role(
        **context.to_dict(),
        target_member_snowflake=member.id,
    )
    server_mute = await database_factory.select(
        singular=True, guild_snowflake=context.guild.id, member_snowflake=member.id
    )
    if not server_mute:
        server_mute = MODEL(
            guild_snowflake=context.guild.id,
            member_snowflake=member.id,
            reason=reason,
        )
        await database_factory.create(server_mute)
        action = "muted"
        should_be_muted = True
    else:
        await database_factory.delete(
            guild_snowflake=context.guild.id, member_snowflake=member.id
        )
        action = "unmuted"
        should_be_muted = False

    if member.voice and member.voice.channel:
        await member.edit(mute=should_be_muted)
    return f"Successfully server {action} {member.mention} in {context.guild.name}."


async def enforce(after, member) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    server_mute = await database_factory.select(
        member_snowflake=member.id, singular=True
    )
    if server_mute:
        if member.guild.id == server_mute.guild_snowflake:
            if not after.mute:
                try:
                    await member.edit(mute=True, reason="Server mute is active.")
                except discord.Forbidden as e:
                    bot.logger.warning(
                        f"No permission to "
                        f"edit mute for {member.display_name}. {str(e).capitalize()}"
                    )
                except discord.HTTPException as e:
                    bot.logger.warning(
                        f"Failed to edit mute for "
                        f"{member.display_name}: "
                        f"{str(e).capitalize()}"
                    )


async def is_server_muted(channel, member) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    server_mute = await database_factory.select(
        channel_snowflake=channel.id, member_snowflake=member.id, singular=False
    )
    if server_mute:
        return True
    return False
