"""!/bin/python3
list_server_mutes.py The purpose of this program is list server mutes.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = VoiceMute


async def build_dictionary(
    guild_snowflake: int, obj
) -> dict[int, dict[str, dict[int, dict[str, dict]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    target = "server"
    server_mutes = []
    dictionary: dict[int, dict[str, dict[int, dict[str, dict]]]] = {}
    if isinstance(obj, discord.Guild):
        server_mutes = await database_factory.select(
            guild_snowflake=obj.id, target=target, singular=False
        )
        guild_snowflake = obj.id
    elif isinstance(obj, discord.Member):
        server_mutes = await database_factory.select(
            guild_snowflake=guild_snowflake,
            member_snowflake=obj.id,
            target=target,
            singular=False,
        )
    else:
        server_mutes = await database_factory.select(
            guild_snowflake=guild_snowflake, target=target, singular=False
        )
    if server_mutes:
        for server_mute in server_mutes:
            dictionary.setdefault(guild_snowflake, {"members": {}})
            dictionary[guild_snowflake]["members"].setdefault(
                server_mute.member_snowflake,
                {"server_mutes": {"reason": server_mute.reason}},
            )
    return dictionary


async def build_pages(guild_snowflake: int, obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return "No active server-mutes found."

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
            return "No active server-mutes found."

    title = f"{emojis.get_random_emoji()} Server Mutes for {obj_name}"

    dictionary = await build_dictionary(guild_snowflake=guild_snowflake, obj=obj)
    processed_dictionary: list_service.ServerMuteDictionary = (
        await list_service.process_dictionary(
            cls=list_service.ServerMuteDictionary, dictionary=dictionary
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        smute_n = 0
        field_count = 0
        lines = []
        thumbnail = False
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for member_snowflake, dictionary in guild_data.get("members", {}).items():
            member = guild.get_member(member_snowflake)
            if member:
                if not isinstance(obj, discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                    lines.append(
                        f"**Reason:** {dictionary.get("server_mutes", {}).get('reason', None)}"
                    )
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
                    lines.append(f"**User:** {display_name} {member_snowflake}")
                else:
                    continue
            smute_n += 1
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
        embed.description = f"**{original_description} ({smute_n})**"
        pages.append(embed)
    if not pages:
        return "No server mutes found."
    return pages
