"""server_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the server mute moderation.

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
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_dictionary,
    flush_page,
)
from vyrtuous.utils.check import (
    has_equal_or_lower_role,
)
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.db.infractions.server_mute import ServerMute


class ServerMuteService:

    lines, pages = [], []

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        server_mutes = await ServerMute.select(singular=False, **where_kwargs)
        for server_mute in server_mutes:
            dictionary.setdefault(server_mute.guild_snowflake, {"members": {}})
            dictionary[server_mute.guild_snowflake]["members"].setdefault(
                server_mute.member_snowflake, {"server_mutes": {}}
            )
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                ServerMute.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                ServerMute.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} {ServerMute.PLURAL} {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await ServerMute.build_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, server_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    ServerMute.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(ServerMute.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, ServerMute.pages, title, guild.name)
                    ServerMute.lines = []
                    field_count = 0
            if ServerMute.lines:
                embed.add_field(
                    name="Information", value="\n".join(ServerMute.lines), inline=False
                )
            ServerMute.pages.append(embed)
        return ServerMute.pages

    @classmethod
    async def toggle_server_mute(cls, member_dict, reason, snowflake_kwargs):
        bot = DiscordBot.get_instance()
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        guild = bot.get_guild(guild_snowflake)
        await has_equal_or_lower_role(
            snowflake_kwargs=snowflake_kwargs,
            member_snowflake=member_dict.get("id", None),
        )
        where_kwargs = member_dict.get("columns", None)

        server_mute = await ServerMute.select(singular=True, **where_kwargs)
        if not server_mute:
            server_mute = ServerMute(
                **where_kwargs,
                reason=reason,
            )
            await server_mute.create()
            action = "muted"
            should_be_muted = True
        else:
            await ServerMute.delete(**where_kwargs)
            action = "unmuted"
            should_be_muted = False

        if (
            member_dict.get("object", None).voice
            and member_dict.get("object", None).voice.channel
        ):
            await member_dict.get("object", None).edit(mute=should_be_muted)
        return f"Successfully server {action} {member_dict.get("mention", None)} in {guild.name}."
