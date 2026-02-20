"""!/bin/python3
server_mute_service.py The purpose of this program is to extend AliasService to service the server mute infraction.

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

from vyrtuous.server_mute.server_mute import ServerMute


@dataclass
class ServerMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, dict]]]] = field(default_factory=dict)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class ServerMuteService:
    __CHUNK_SIZE = 7
    MODEL = ServerMute

    def __init__(
        self,
        *,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
        moderator_service=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__moderator_service = moderator_service

    async def build_dictionary(self, where_kwargs):
        dictionary = {}
        server_mutes = await self.__database_factory.select(
            singular=False, **where_kwargs
        )
        for server_mute in server_mutes:
            dictionary.setdefault(server_mute.guild_snowflake, {"members": {}})
            dictionary[server_mute.guild_snowflake]["members"].setdefault(
                server_mute.member_snowflake, {"server_mutes": {}}
            )
        return dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Server Mutes {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_dictionary(where_kwargs=where_kwargs)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=ServerMuteDictionary, dictionary=dictionary
        )
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)

        smute_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            field_count = 0
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, server_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                smute_n += 1
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
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({smute_n})**"
        return pages

    async def toggle_server_mute(self, context, member_dict, reason):
        guild_snowflake = context.guild.id
        guild = self.__bot.get_guild(guild_snowflake)
        await self.__moderator_service.has_equal_or_lower_role(
            **context.to_dict(),
            target_member_snowflake=member_dict.get("id", None),
        )
        server_mute = await self.__database_factory.select(
            singular=True, **member_dict.get("columns", None)
        )
        if not server_mute:
            server_mute = self.MODEL(
                **member_dict.get("columns", None),
                reason=reason,
            )
            await self.__database_factory.create(server_mute)
            action = "muted"
            should_be_muted = True
        else:
            await self.__database_factory.delete(**member_dict.get("columns", None))
            action = "unmuted"
            should_be_muted = False

        if (
            member_dict.get("object", None).voice
            and member_dict.get("object", None).voice.channel
        ):
            await member_dict.get("object", None).edit(mute=should_be_muted)
        return f"Successfully server {action} {member_dict.get('mention', None)} in {guild.name}."

    async def enforce(self, after, member):
        server_mute = await self.__database_factory.select(
            member_snowflake=member.id, singular=True
        )
        if server_mute:
            if member.guild.id == server_mute.guild_snowflake:
                if not after.mute:
                    try:
                        await member.edit(mute=True, reason="Server mute is active.")
                    except discord.Forbidden as e:
                        self.__bot.logger.warning(
                            f"No permission to "
                            f"edit mute for {member.display_name}. {str(e).capitalize()}"
                        )
                    except discord.HTTPException as e:
                        self.__bot.logger.warning(
                            f"Failed to edit mute for "
                            f"{member.display_name}: "
                            f"{str(e).capitalize()}"
                        )
            return False
        return True

    async def is_server_muted(self, channel, member):
        server_mute = await self.__database_factory.select(
            channel_snowflake=channel.id, member_snowflake=member.id
        )
        if server_mute:
            return True
        return False
