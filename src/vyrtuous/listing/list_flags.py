"""!/bin/python3
list_flags.py The purpose of this program is list flags.

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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.flag import Flag
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = Flag


async def build_dictionary(
    guild_snowflake: int,
    obj,
) -> dict[int, dict[str, dict[int, dict[str, dict[int, dict[str, str]]]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    flags = []
    dictionary: dict[
        int, dict[str, dict[int, dict[str, dict[int, dict[str, str]]]]]
    ] = {}
    if isinstance(obj, discord.Guild):
        flags = await database_factory.select(guild_snowflake=obj.id, singular=False)
        guild_snowflake = obj.id
    elif isinstance(obj, discord.abc.GuildChannel):
        flags = await database_factory.select(
            channel_snowflake=obj.id, guild_snowflake=guild_snowflake, singular=False
        )
    elif isinstance(obj, discord.Member):
        flags = await database_factory.select(
            member_snowflake=obj.id, guild_snowflake=guild_snowflake, singular=False
        )
    elif isinstance(obj, int):
        flags = await database_factory.select(
            member_snowflake=obj, guild_snowflake=guild_snowflake, singular=False
        )
    else:
        flags = await database_factory.select(
            guild_snowflake=guild_snowflake, singular=False
        )
    if flags:
        for flag in flags:
            dictionary.setdefault(guild_snowflake, {"members": {}})
            dictionary[guild_snowflake]["members"].setdefault(
                flag.member_snowflake, {"flags": {}}
            )
            dictionary[guild_snowflake]["members"][flag.member_snowflake][
                "flags"
            ].setdefault(flag.channel_snowflake, {})
            dictionary[guild_snowflake]["members"][flag.member_snowflake]["flags"][
                flag.channel_snowflake
            ].update(
                {
                    "reason": flag.reason,
                }
            )
    return dictionary


async def build_pages(guild_snowflake: int, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return "No active flags found."
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
    title = f"{emojis.get_random_emoji()} Flags for {obj_name}"

    dictionary = await build_dictionary(guild_snowflake=guild_snowflake, obj=obj)
    processed_dictionary: list_service.FlagDictionary = (
        await list_service.process_dictionary(
            cls=list_service.FlagDictionary, dictionary=dictionary
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        flag_n = 0
        field_count = 0
        lines = []
        thumbnail = False
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for member_snowflake, flag_dictionary in guild_data.get("members", {}).items():
            member = guild.get_member(member_snowflake)
            if member:
                if not isinstance(obj, discord.Member):
                    lines.append(f"**User:** {member.mention}")
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
            flag_n += 1
            for channel_snowflake, channel_dictionary in flag_dictionary.get(
                "flags", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    continue
                if not isinstance(obj, discord.abc.GuildChannel):
                    lines.append(f"**Channel:** {channel.mention}")
                if isinstance(obj, discord.Member):
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
        if lines:
            embed.add_field(name="Information", value="\n".join(lines), inline=False)
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({flag_n})**"
        pages.append(embed)
    if not pages:
        return "No flags found."
    return pages
