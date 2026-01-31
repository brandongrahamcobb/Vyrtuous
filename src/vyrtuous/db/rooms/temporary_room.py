"""temporary_rooms.py A utility module for managing temporary rooms in the Vyrtuous Discord bot.

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

from datetime import datetime
from typing import Optional

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_dictionary,
    flush_page,
)

from vyrtuous.db.actions.ban import Ban
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.db.roles.coordinator import Coordinator
from vyrtuous.db.actions.flag import Flag
from vyrtuous.db.roles.moderator import Moderator
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.db.rooms.stage import Stage


class TemporaryRoom(DatabaseFactory):

    ACT = "temp"

    CATEGORY = "temp"
    PLURAL = "Temporary Rooms"
    SCOPES = ["channels"]
    SINGULAR = "Temporary Rooms"
    UNDO = "temp"

    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "room_name",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "updated_at",
    ]

    TABLE_NAME = "temporary_rooms"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        room_name: str,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_mention = f"<#{channel_snowflake}>"
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.is_temp_room: Optional[bool] = True
        self.room_name = room_name
        self.updated_at = updated_at


    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        chunk_size, field_count, pages = 7, 0, []
        dictionary = {}
        aliases = await Alias.select(**where_kwargs)
        temporary_rooms = await TemporaryRoom.select(**where_kwargs)
        for temporary_room in temporary_rooms:
            dictionary.setdefault(
                temporary_room.guild_snowflake, {"channels": {}}
            )
            dictionary[temporary_room.guild_snowflake]["channels"].setdefault(
                temporary_room.channel_snowflake, {}
            )
            if aliases:
                for alias in aliases:
                    if (
                        alias.guild_snowflake == temporary_room.guild_snowflake
                        and alias.channel_snowflake == temporary_room.channel_snowflake
                    ):
                        dictionary[temporary_room.guild_snowflake]["channels"][
                            temporary_room.channel_snowflake
                        ].setdefault(alias.category, [])
                        dictionary[temporary_room.guild_snowflake]["channels"][
                            temporary_room.channel_snowflake
                        ][alias.category].append(alias.alias_name)
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
        return cleaned_dictionary, pages
        
    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines = 7, 0, []
        title = f"{get_random_emoji()} {TemporaryRoom.PLURAL}"

        where_kwargs = object_dict.get("columns", None)

        dictionary, pages = await TemporaryRoom.build_clean_dictionary(is_at_home=is_at_home, where_kwargs=where_kwargs)

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for category, alias_names in channel_data.items():
                    lines.append(f"{category}")
                    field_count += 1
                    for name in alias_names:
                        lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= chunk_size:
                            embed.add_field(
                                name="Information",
                                value="\n".join(lines),
                                inline=False,
                            )
                            embed, field_count = flush_page(
                                embed, pages, title, guild.name
                            )
                            lines = []
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)
        return pages

    @classmethod
    async def migrate_temporary_room(cls, channel_dict, old_name, snowflake_kwargs):
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        old_room = await TemporaryRoom.select(
            guild_snowflake=guild_snowflake, room_name=old_name, singular=True
        )
        set_kwargs = {"channel_snowflake": channel_dict.get("id", None)}
        temp_where_kwargs = {
            "channel_snowflake": old_room.channel_snowflake,
            "guild_snowflake": guild_snowflake,
            "room_name": channel_dict.get("name", None),
        }
        where_kwargs = {
            "channel_snowflake": old_room.channel_snowflake,
            "guild_snowflake": guild_snowflake,
        }
        kwargs = {
            "set_kwargs": set_kwargs,
            "where_kwargs": where_kwargs,
        }
        await TemporaryRoom.update(
            set_kwargs=set_kwargs,
            where_kwargs=temp_where_kwargs,
        )
        await Alias.update(**kwargs)
        await Ban.update(**kwargs)
        await Cap.update(**kwargs)
        await Coordinator.update(**kwargs)
        await Flag.update(**kwargs)
        await Moderator.update(**kwargs)
        await Stage.update(**kwargs)
        await TextMute.update(**kwargs)
        await VoiceMute.update(**kwargs)
        await Vegan.update(**kwargs)
        return f"Temporary room `{old_name}` migrated to {channel_dict.get("mention", None)}."

    @classmethod
    async def toggle_temporary_room(cls, channel_dict):
        action = None
        kwargs = {}
        kwargs.update(channel_dict.get("columns", None))
        temporary_room = await TemporaryRoom.select(**kwargs, singular=True)
        if temporary_room:
            await TemporaryRoom.delete(**kwargs)
            action = "removed"
        else:
            temporary_room = TemporaryRoom(
                **kwargs,
                room_name=channel_dict.get("name", None),
            )
            await temporary_room.create()
            action = "created"
        return f"Temporary room {action} in {channel_dict.get("mention", None)}."
