"""!/bin/python3
alias_service.py The purpose of this program is to extend Service to service aliases.

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
from vyrtuous.db.alias import Alias
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.listing import list_service
from vyrtuous.utils.messaging import emojis

MODEL = Alias


async def build_dictionary(obj) -> dict:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    aliases = []
    dictionary: dict[
        int, dict[str, dict[int, dict[str, dict[str, dict[str, list | str]]]]]
    ] = {}
    if isinstance(obj, discord.Guild):
        aliases = await database_factory.select(guild_snowflake=obj.id, singular=False)
    elif isinstance(obj, discord.abc.GuildChannel):
        aliases = await database_factory.select(
            channel_snowflake=obj.id, singular=False
        )
    else:
        aliases = await database_factory.select(singular=False)
    if aliases:
        for alias in aliases:
            dictionary.setdefault(alias.guild_snowflake, {"channels": {}})
            dictionary[alias.guild_snowflake]["channels"].setdefault(
                alias.channel_snowflake, {"aliases": {}}
            )
            dictionary[alias.guild_snowflake]["channels"][alias.channel_snowflake][
                "aliases"
            ].setdefault(alias.category, {})[alias.alias_name] = []
            if alias.category == "role":
                guild = bot.get_guild(alias.guild_snowflake)
                if guild:
                    role = guild.get_role(alias.role_snowflake)
                    if not role:
                        continue
                    dictionary[alias.guild_snowflake]["channels"][
                        alias.channel_snowflake
                    ]["aliases"][alias.category][alias.alias_name] = role.mention
    return dictionary


async def build_pages(obj) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Command Aliases in {obj_name}"

    dictionary = await build_dictionary(obj=obj)
    processed_dictionary: list_service.AliasDictionary = (
        await list_service.process_dictionary(
            cls=list_service.AliasDictionary, dictionary=dictionary
        )
    )

    for guild_snowflake, guild_data in processed_dictionary.data.items():
        alias_n = 0
        field_count = 0
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            continue
        embed = discord.Embed(
            title=title, description=guild.name, color=discord.Color.blue()
        )
        for channel_snowflake, channel_dictionary in guild_data.get(
            "channels", {}
        ).items():
            channel = guild.get_channel(channel_snowflake)
            if channel is None:
                continue
            for category, alias_data in channel_dictionary.get("aliases", {}).items():
                lines.append(f"{category}")
                for name, role_mention in alias_data.items():
                    if category == "role":
                        lines.append(f"  ↳ `{name}` -> {role_mention}")
                    else:
                        lines.append(f"  ↳ `{name}`")
                    alias_n += 1
                field_count += 1
            if field_count >= list_service.CHUNK_SIZE:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
                embed = list_service.flush_page(embed, pages, title, guild.name)
                lines = []
                field_count = 0
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({alias_n})**"
        pages.append(embed)
    if not pages:
        return "No aliases found."
    return pages
