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

from copy import copy
from dataclasses import dataclass, field
from typing import Dict, List

import discord

from vyrtuous.temporary_room.temporary_room import TemporaryRoom


@dataclass
class TemporaryRoomDictionary:
    data: Dict[int, Dict[str, Dict[int, bool]]] = field(default_factory=dict)
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


class TemporaryRoomService:
    __CHUNK_SIZE = 12
    MODEL = TemporaryRoom

    def __init__(
        self,
        *,
        alias_service=None,
        bot=None,
        cap_service=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
        moderator_service=None,
        stage_service=None,
        coordinator_service=None,
        voice_mute_service=None,
        ban_service=None,
        flag_service=None,
        vegan_service=None,
        text_mute_service=None,
    ):
        self.__alias_service = alias_service
        self.__bot = bot
        self.__cap_service = cap_service
        self.__database_factory = copy(database_factory)
        self.__dictionary_service = dictionary_service
        self.__database_factory.model = self.MODEL
        self.__emoji = emoji
        self.__moderator_service = moderator_service
        self.__stage_service = stage_service
        self.__coordinator_service = coordinator_service
        self.__voice_mute_service = voice_mute_service
        self.__ban_service = ban_service
        self.__flag_service = flag_service
        self.__vegan_service = vegan_service
        self.__text_mute_service = text_mute_service
        self.deleted_rooms = {}

    async def build_dictionary(self, obj):
        temporary_rooms = []
        dictionary = {}
        if isinstance(obj, discord.Guild):
            temporary_rooms = await self.__database_factory.select(
                guild_snowflake=obj.id
            )
        elif isinstance(obj, discord.abc.GuildChannel):
            temporary_rooms = await self.__database_factory.select(
                channel_snowflake=obj.id
            )
        else:
            temporary_rooms = await self.__database_factory.select()
        if temporary_rooms:
            for temporary_room in temporary_rooms:
                dictionary.setdefault(temporary_room.guild_snowflake, {"channels": {}})
                dictionary[temporary_room.guild_snowflake]["channels"][
                    temporary_room.channel_snowflake
                ] = True
        return dictionary

    async def build_pages(self, is_at_home, obj):
        lines, pages = [], []

        obj_name = "All Servers"
        if obj:
            obj_name = obj.name
        title = f"{self.__emoji.get_random_emoji()} Temporary Rooms for {obj_name}"

        dictionary = await self.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=TemporaryRoomDictionary, dictionary=dictionary
        )

        for guild_snowflake, guild_data in processed_dictionary.data.items():
            temp_n = 0
            field_count = 0
            lines = []
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for category, alias_names in channel_data.items():
                    lines.append(f"{category}")
                    temp_n += 1
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= self.__CHUNK_SIZE:
                            embed.add_field(
                                name="Information",
                                value="\n".join(lines),
                                inline=False,
                            )
                            embed = self.__dictionary_service.flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                            field_count = 0
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            original_description = embed.description or ""
            embed.description = f"**{original_description}** **({temp_n})**"
            pages.append(embed)
        if is_at_home:
            pages.extend(processed_dictionary.skipped_channels)
            pages.extend(processed_dictionary.skipped_guilds)
        return pages

    async def migrate_temporary_room(self, channel, old_name):
        old_room = await self.__database_factory.select(
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
        await self.__database_factory.update(
            set_kwargs=set_kwargs,
            where_kwargs=temp_where_kwargs,
        )
        await self.__alias_service.migrate(kwargs)
        await self.__ban_service.migrate(kwargs)
        await self.__cap_service.migrate(kwargs)
        await self.__coordinator_service.migrate(kwargs)
        await self.__flag_service.migrate(kwargs)
        await self.__moderator_service.migrate(kwargs)
        await self.__stage_service.migrate(kwargs)
        await self.__text_mute_service.migrate(kwargs)
        await self.__voice_mute_service.migrate(kwargs)
        await self.__vegan_service.migrate(kwargs)
        return f"Temporary room `{old_name}` migrated to {channel.mention}."

    async def toggle_temporary_room(self, channel):
        temporary_room = await self.__database_factory.select(
            channel_snowflake=channel.id, singular=True
        )
        if temporary_room:
            await self.__database_factory.delete(channel_snowflake=channel.id)
            action = "removed"
        else:
            temporary_room = self.MODEL(
                channel_snowflake=int(channel.id),
                guild_snowflake=int(channel.guild.id),
                room_name=str(channel.name),
            )
            await self.__database_factory.create(temporary_room)
            action = "created"
        return f"Temporary room {action} in {channel.mention}."

    async def add_deleted_room(self, channel):
        room = await self.__database_factory.select(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            singular=True,
        )
        if room:
            self.deleted_rooms[channel.name] = room

    async def rename_room(self, before, after):
        set_kwargs = {"room_name": after.id}
        where_kwargs = {"room_name": before.id}
        await self.__database_factory.update(
            set_kwargs=set_kwargs, where_kwargs=where_kwargs
        )
