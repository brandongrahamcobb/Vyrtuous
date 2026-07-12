"""!/bin/python3
coordinator_service.py The purpose of this program is to extend Service to service the coordinator class.

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
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.coordinator import Coordinator
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = Coordinator


async def build_dictionary(
    obj,
) -> dict[int, dict[str, dict[int, dict[str, dict[int, dict[str, str]]]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    coordinators = []
    dictionary: dict[
        int, dict[str, dict[int, dict[str, dict[int, dict[str, str]]]]]
    ] = {}
    if isinstance(obj, discord.Guild):
        coordinators = await database_factory.select(
            guild_snowflake=obj.id, singular=False
        )
    elif isinstance(obj, discord.abc.GuildChannel):
        coordinators = await database_factory.select(
            channel_snowflake=obj.id, singular=False
        )
    elif isinstance(obj, discord.Member):
        coordinators = await database_factory.select(
            member_snowflake=obj.id, singular=False
        )
    else:
        coordinators = await database_factory.select(singular=False)
    if coordinators:
        for coordinator in coordinators:
            dictionary.setdefault(coordinator.guild_snowflake, {"members": {}})
            dictionary[coordinator.guild_snowflake]["members"].setdefault(
                coordinator.member_snowflake, {"coordinators": {}}
            )
            dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"].setdefault(coordinator.channel_snowflake, {})
            dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"][coordinator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )
    return dictionary


async def build_pages(is_at_home: bool, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = "All Servers"
    if obj is not None and not isinstance(obj, (int, str)):
        obj_name = obj.name
    elif isinstance(obj, int):
        simplified_member = bot.registry.get(MemberState).active.get(obj, None)
        if simplified_member:
            obj_name = simplified_member[0]
        else:
            return "No coordinators found."
    title = f"{emojis.get_random_emoji()} Coordinators for {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary: list_service.CoordinatorDictionary = (
        await list_service.process_dictionary(
            cls=list_service.CoordinatorDictionary, dictionary=dictionary
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        coord_n = 0
        field_count = 0
        lines = []
        thumbnail = False
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for member_snowflake, coordinator_dictionary in guild_data.get(
            "members", {}
        ).items():
            member = guild.get_member(member_snowflake)
            if member:
                if not isinstance(obj, discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
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
            coord_n += 1
            for channel_snowflake, _ in coordinator_dictionary.get(
                "coordinators", {}
            ).items():
                if not isinstance(obj, discord.abc.GuildChannel):
                    channel = guild.get_channel(channel_snowflake)
                    if channel is None:
                        continue
                    lines.append(f"**Channel:** {channel.mention}")
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
        embed.description = f"**{original_description} ({coord_n})**"
        pages.append(embed)
    if is_at_home:
        pages.extend(processed_dictionary.skipped_guilds)
        pages.extend(processed_dictionary.skipped_members)
    if not pages:
        return "No coordinators found."
    return pages
