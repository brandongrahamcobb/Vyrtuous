"""!/bin/python3
list_text_mutes.py The purpose of this program is to list text-mutes.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from typing import Any

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.text_mute import TextMute
from vyrtuous.listing import list_service
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.messaging import emojis

MODEL = TextMute


async def build_dictionary(
    guild_snowflake: int,
    obj,
) -> dict[int, dict[str, dict[int, dict[str, dict[int, dict[str, Any]]]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    text_mutes = []
    dictionary: dict[
        int, dict[str, dict[int, dict[str, dict[int, dict[str, Any]]]]]
    ] = {}
    if isinstance(obj, discord.Guild):
        text_mutes = await database_factory.select(
            guild_snowflake=obj.id, singular=False
        )
        guild_snowflake = obj.id
    elif isinstance(obj, discord.abc.GuildChannel):
        text_mutes = await database_factory.select(
            channel_snowflake=obj.id, guild_snowflake=guild_snowflake, singular=False
        )
    elif isinstance(obj, discord.Member):
        text_mutes = await database_factory.select(
            guild_snowflake=guild_snowflake, member_snowflake=obj.id, singular=False
        )
    elif isinstance(obj, int):
        text_mutes = await database_factory.select(
            guild_snowflake=guild_snowflake, member_snowflake=obj, singular=False
        )
    else:
        text_mutes = await database_factory.select(
            guild_snowflake=guild_snowflake, singular=False
        )
    if text_mutes:
        for text_mute in text_mutes:
            dictionary.setdefault(guild_snowflake, {"members": {}})
            dictionary[guild_snowflake]["members"].setdefault(
                text_mute.member_snowflake, {"text_mutes": {}}
            )
            dictionary[guild_snowflake]["members"][text_mute.member_snowflake][
                "text_mutes"
            ].setdefault(text_mute.channel_snowflake, {})
            dictionary[guild_snowflake]["members"][text_mute.member_snowflake][
                "text_mutes"
            ][text_mute.channel_snowflake].update(
                {"reason": text_mute.reason, "expires_in": text_mute.expires_in}
            )
    return dictionary


async def build_pages(guild_snowflake: int, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return "No active text-mutes found."
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = guild.name
    if not isinstance(obj, int):
        obj_name = obj.name
    else:
        simplified_member = bot.registry.get(MemberState).active.get(obj, None)
        if simplified_member:
            obj_name = simplified_member[0]
        else:
            return "This command must target a valid member."

    title = f"{emojis.get_random_emoji()} Text Mutes for {obj_name}"

    dictionary = await build_dictionary(guild_snowflake=guild_snowflake, obj=obj)
    processed_dictionary: list_service.TextMuteDictionary = (
        await list_service.process_dictionary(
            cls=list_service.TextMuteDictionary, dictionary=dictionary
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        duration_builder: DurationBuilder = DurationBuilder()
        tmute_n = 0
        field_count = 0
        lines = []
        thumbnail = False
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for member_snowflake, text_mute_dictionary in guild_data.get(
            "members", {}
        ).items():
            member = guild.get_member(member_snowflake)
            if member:
                if not isinstance(obj, discord.Member):
                    lines.append(f"**User:** {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(url=obj.display_avatar.url)
                    thumbnail = True
            else:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake, None
                )
                if simplified_member:
                    display_name = simplified_member[0]
                    lines.append(f"**User:** {display_name} ({member_snowflake})")
                else:
                    continue
            tmute_n += 1
            for channel_snowflake, channel_dictionary in text_mute_dictionary.get(
                "text_mutes", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    continue
                if not isinstance(obj, discord.abc.GuildChannel):
                    lines.append(f"**Channel:** {channel.mention}")
                if isinstance(obj, discord.Member):
                    lines.append(
                        f"**Expires:** {duration_builder.from_timestamp(channel_dictionary['expires_in']).to_unix_ts()}"
                    )
                    lines.append(f"**Reason:** {channel_dictionary['reason']}")
                field_count += 1
                if field_count >= list_service.CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = list_service.flush_page(embed, pages, title, guild.name)
                    lines = []
                    field_count = 0
        if lines:
            embed.add_field(
                name="Information",
                value="\n".join(lines),
                inline=False,
            )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({tmute_n})**"
        pages.append(embed)
    if not pages:
        return "No text mutes found."
    return pages
