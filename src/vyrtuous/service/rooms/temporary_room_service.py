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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.infractions.ban import Ban
from vyrtuous.db.infractions.flag import Flag
from vyrtuous.db.infractions.text_mute import TextMute
from vyrtuous.db.infractions.voice_mute import VoiceMute
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.db.roles.coordinator import Coordinator
from vyrtuous.db.roles.moderator import Moderator
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.db.rooms.stage import Stage
from vyrtuous.db.rooms.temporary_room import TemporaryRoom
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_channels,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji


class TemporaryRoomService:

    lines, pages = [], []

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        aliases = await Alias.select(singular=False, **where_kwargs)
        temporary_rooms = await TemporaryRoom.select(singular=False, **where_kwargs)
        for temporary_room in temporary_rooms:
            dictionary.setdefault(temporary_room.guild_snowflake, {"channels": {}})
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
                TemporaryRoom.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
            if skipped_guilds:
                TemporaryRoom.pages.extend(
                    pages=generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Temporary Rooms"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await TemporaryRoomService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

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
                TemporaryRoom.lines.append(f"Channel: {channel.mention}")
                field_count += 1
                for category, alias_names in channel_data.items():
                    TemporaryRoom.lines.append(f"{category}")
                    field_count += 1
                    for name in alias_names:
                        TemporaryRoom.lines.append(f"  ↳ {name}")
                        field_count += 1
                        if field_count >= CHUNK_SIZE:
                            embed.add_field(
                                name="Information",
                                value="\n".join(TemporaryRoom.lines),
                                inline=False,
                            )
                            embed = flush_page(
                                embed, TemporaryRoom.pages, title, guild.name
                            )
                            TemporaryRoom.lines = []
                            field_count = 0
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(TemporaryRoom.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, TemporaryRoom.pages, title, guild.name)
                    TemporaryRoom.lines = []
                    field_count = 0
            if TemporaryRoom.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(TemporaryRoom.lines),
                    inline=False,
                )
            TemporaryRoom.pages.append(embed)
        return TemporaryRoom.pages

    @classmethod
    async def migrate_temporary_room(cls, channel_dict, old_name, snowflake_kwargs):
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        old_room = await TemporaryRoom.select(
            guild_snowflake=int(guild_snowflake), room_name=old_name, singular=True
        )
        if not old_room:
            return f"No temporary room found with the name {old_name}."
        set_kwargs = {"channel_snowflake": channel_dict.get("id", None)}
        temp_where_kwargs = {
            "channel_snowflake": old_room.channel_snowflake,
            "guild_snowflake": int(guild_snowflake),
            "room_name": channel_dict.get("name", None),
        }
        where_kwargs = {
            "channel_snowflake": old_room.channel_snowflake,
            "guild_snowflake": int(guild_snowflake),
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
