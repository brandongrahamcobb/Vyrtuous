"""!/bin/python3
temporary_rooms_service.py The purpose of this program is to extend Service to service the temporary room class.

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

from vyrtuous.aliases import alias_service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.temporary_room import TemporaryRoom
from vyrtuous.utils.moderation import (ban_service, flag_service,
                                       text_mute_service, voice_mute_service)
from vyrtuous.utils.rooms import automute_room_service, cap_service
from vyrtuous.utils.users import coordinator_service, moderator_service
from vyrtuous.vegan import vegan_service

MODEL = TemporaryRoom


async def migrate_temporary_room(channel, old_name):
    database_factory = DatabaseFactory(MODEL)
    old_room = await database_factory.select(
        guild_snowflake=int(channel.guild.id),
        room_name=str(old_name),
        singular=True,
    )
    if not old_room:
        return f"No temporary room found with the name {old_name}."
    set_kwargs = {"channel_snowflake": channel.id}
    temp_where_kwargs = {
        "channel_snowflake": int(old_room.channel_snowflake),
        "guild_snowflake": int(channel.guild.id),
        "room_name": str(channel.name),
    }
    where_kwargs = {
        "channel_snowflake": int(old_room.channel_snowflake),
        "guild_snowflake": int(channel.guild.id),
    }
    kwargs = {
        "set_kwargs": set_kwargs,
        "where_kwargs": where_kwargs,
    }
    await database_factory.update(
        set_kwargs=set_kwargs,
        where_kwargs=temp_where_kwargs,
    )
    await alias_service.migrate(kwargs)
    await ban_service.migrate(kwargs)
    await cap_service.migrate(kwargs)
    await coordinator_service.migrate(kwargs)
    await flag_service.migrate(kwargs)
    await moderator_service.migrate(kwargs)
    await automute_room_service.migrate(kwargs)
    await text_mute_service.migrate(kwargs)
    await voice_mute_service.migrate(kwargs)
    await vegan_service.migrate(kwargs)
    return f"Temporary room `{old_name}` migrated to {channel.mention}."


async def toggle_temporary_room(channel):
    database_factory = DatabaseFactory(MODEL)
    temporary_room = await database_factory.select(
        channel_snowflake=channel.id, singular=True
    )
    if temporary_room:
        await database_factory.delete(channel_snowflake=channel.id)
        action = "removed"
    else:
        temporary_room = MODEL(
            channel_snowflake=int(channel.id),
            guild_snowflake=int(channel.guild.id),
            room_name=str(channel.name),
        )
        await database_factory.create(temporary_room)
        action = "created"
    return f"Temporary room {action} in {channel.mention}."


async def add_deleted_room(channel):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    room = await database_factory.select(
        channel_snowflake=channel.id,
        guild_snowflake=channel.guild.id,
        singular=True,
    )
    if room:
        bot.registry.get(ChannelState).deleted[channel.name] = room.channel_snowflake


async def rename_room(before, after):
    database_factory = DatabaseFactory(MODEL)
    set_kwargs = {"room_name": after.id}
    where_kwargs = {"room_name": before.id}
    await database_factory.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
