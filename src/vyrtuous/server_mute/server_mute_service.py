from copy import copy

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

import discord

from vyrtuous.server_mute.server_mute import ServerMute


class ServerMuteService:
    __CHUNK_SIZE = 7
    MODEL = ServerMute

    def __init__(
        self, *, bot=None, database_factory=None, dictionary_service=None, emoji=None
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        dictionary = {}
        server_mutes = await self.__database_factory.select(
            singular=False, **where_kwargs
        )
        for server_mute in server_mutes:
            dictionary.setdefault(server_mute.guild_snowflake, {"members": {}})
            dictionary[server_mute.guild_snowflake]["members"].setdefault(
                server_mute.member_snowflake, {"server_mutes": {}}
            )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        skipped_members = self.__dictionary_service.generate_skipped_members(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        # if is_at_home:
        #     if skipped_guilds:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_set_pages(
        #                 skipped=skipped_guilds,
        #                 title="Skipped Servers",
        #             )
        #         )
        #     if skipped_members:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_dict_pages(
        #                 skipped=skipped_members,
        #                 title="Skipped Members in Server",
        #             )
        #         )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Server Mutes {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )
        embed = discord.Embed(
            title=title, description="Default view", color=discord.Color.blue()
        )
        if not dictionary:
            pages = [embed]

        smute_n = 0
        for guild_snowflake, guild_data in dictionary.items():
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

    async def toggle_server_mute(self, default_kwargs, member_dict, reason):
        updated_kwargs = default_kwargs.copy()
        guild_snowflake = default_kwargs.get("guild_snowflake", None)
        guild = self.__bot.get_guild(guild_snowflake)
        # await PermissionService.has_equal_or_lower_role(
        #     **default_kwargs,
        #     target_member_snowflake=member_dict.get("id", None),
        # )
        updated_kwargs.update(member_dict.get("columns", None))
        del updated_kwargs["channel_snowflake"]
        server_mute = await self.__database_factory.select(
            singular=True, **updated_kwargs
        )
        if not server_mute:
            server_mute = self.MODEL(
                **updated_kwargs,
                reason=reason,
            )
            await self.__database_factory.create(server_mute)
            action = "muted"
            should_be_muted = True
        else:
            await self.__database_factory.delete(**updated_kwargs)
            action = "unmuted"
            should_be_muted = False

        if (
            member_dict.get("object", None).voice
            and member_dict.get("object", None).voice.channel
        ):
            await member_dict.get("object", None).edit(mute=should_be_muted)
        return f"Successfully server {action} {member_dict.get('mention', None)} in {guild.name}."
