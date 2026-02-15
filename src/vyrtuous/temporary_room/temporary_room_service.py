from copy import copy

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

import discord

from vyrtuous.alias.alias import Alias
from vyrtuous.ban.ban import Ban
from vyrtuous.base.record_service import RecordService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cap.cap import Cap
from vyrtuous.coordinator.coordinator import Coordinator
from vyrtuous.flag.flag import Flag
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.moderator.moderator import Moderator
from vyrtuous.stage_room.stage import Stage
from vyrtuous.temporary_room.temporary_room import TemporaryRoom
from vyrtuous.text_mute.text_mute import TextMute
from vyrtuous.vegan.vegan import Vegan
from vyrtuous.voice_mute.voice_mute import VoiceMute


class TemporaryRoomService(RecordService):
    MODEL = TemporaryRoom

    def __init__(
        self,
        *,
        alias_service=None,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
    ):
        self.__alias_service = alias_service
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__dictionary_service = dictionary_service
        self.__database_factory.model = self.MODEL
        self.__emoji = emoji

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        pages = []
        dictionary = {}
        # aliases = await Alias.select(singular=False, **where_kwargs)
        temporary_rooms = await self.__database_factory.select(
            singular=False, **where_kwargs
        )
        for temporary_room in temporary_rooms:
            dictionary.setdefault(temporary_room.guild_snowflake, {"channels": {}})
            dictionary[temporary_room.guild_snowflake]["channels"].setdefault(
                temporary_room.channel_snowflake, {}
            )
            # if aliases:
            #     for alias in aliases:
            #         if (
            #             alias.guild_snowflake == temporary_room.guild_snowflake
            #             and alias.channel_snowflake == temporary_room.channel_snowflake
            #         ):
            #             dictionary[temporary_room.guild_snowflake]["channels"][
            #                 temporary_room.channel_snowflake
            #             ].setdefault(alias.category, [])
            #             dictionary[temporary_room.guild_snowflake]["channels"][
            #                 temporary_room.channel_snowflake
            #             ][alias.category].append(alias.alias_name)
        skipped_channels = self.__dictionary_service.generate_skipped_channels(
            dictionary
        )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_channels:
                pages.extend(
                    self.__dictionary_service.generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
            if skipped_guilds:
                pages.extend(
                    self.__dictionary_service.generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = []
        title = f"{self.__emoji.get_random_emoji()} Temporary Rooms"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        temp_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
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
                        if field_count >= CHUNK_SIZE:
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
            pages.append(embed)
        if pages:
            pages[0].description = f"**({temp_n})**"
        return pages

    async def migrate_temporary_room(self, channel_dict, default_kwargs, old_name):
        guild_snowflake = default_kwargs.get("guild_snowflake", None)
        # old_room = await self.__database_factory.select(
        #     guild_snowflake=int(guild_snowflake), room_name=old_name, singular=True
        # )
        # if not old_room:
        #     return f"No temporary room found with the name {old_name}."
        # set_kwargs = {"channel_snowflake": channel_dict.get("id", None)}
        # temp_where_kwargs = {
        #     "channel_snowflake": old_room.channel_snowflake,
        #     "guild_snowflake": int(guild_snowflake),
        #     "room_name": channel_dict.get("name", None),
        # }
        # where_kwargs = {
        #     "channel_snowflake": old_room.channel_snowflake,
        #     "guild_snowflake": int(guild_snowflake),
        # }
        # kwargs = {
        #     "set_kwargs": set_kwargs,
        #     "where_kwargs": where_kwargs,
        # }
        # await self.__database_factory.update(
        #     set_kwargs=set_kwargs,
        #     where_kwargs=temp_where_kwargs,
        # )
        # await Alias.update(**kwargs)
        # await Ban.update(**kwargs)
        # await Cap.update(**kwargs)
        # await Coordinator.update(**kwargs)
        # await Flag.update(**kwargs)
        # await Moderator.update(**kwargs)
        # await Stage.update(**kwargs)
        # await TextMute.update(**kwargs)
        # await VoiceMute.update(**kwargs)
        # await Vegan.update(**kwargs)
        # return f"Temporary room `{old_name}` migrated to {channel_dict.get('mention', None)}."

    async def toggle_temporary_room(self, channel_dict):
        kwargs = {}
        kwargs.update(channel_dict.get("columns", None))
        temporary_room = await self.__database_factory.select(**kwargs, singular=True)
        if temporary_room:
            await self.__database_factory.delete(**kwargs)
            action = "removed"
        else:
            temporary_room = TemporaryRoom(
                **kwargs,
                room_name=channel_dict.get("name", None),
            )
            await self.__database_factory.create(temporary_room)
            action = "created"
        return f"Temporary room {action} in {channel_dict.get('mention', None)}."
