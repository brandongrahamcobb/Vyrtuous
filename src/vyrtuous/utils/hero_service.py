"""!/bin/python3
permission_service.py The purpose of this program is to provide the service for deciding whether a member has sufficient permissions.

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

from typing import Dict, Tuple

import discord


class HeroService:
    _invincible_members: Dict[Tuple[int, int], bool] = {}
    __CHUNK_SIZE = 7
    state = False

    def __init__(
        self,
        *,
        ban_service=None,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        emoji=None,
        flag_service=None,
        text_mute_service=None,
        voice_mute_service=None,
    ):
        self.__bot = bot
        self.__database_factory = database_factory
        self.__dictionary_service = dictionary_service
        self.__emoji = emoji
        self.__infractions = [
            ban_service.MODEL,
            flag_service.MODEL,
            text_mute_service.MODEL,
            voice_mute_service.MODEL,
        ]

    async def unrestrict(self, guild_snowflake, member_snowflake):
        guild = self.__bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "member_snowflake": member_snowflake,
        }
        for infraction in self.__infractions:
            objects = self.__database_factory.model = infraction
            if objects:
                match infraction.identifier:
                    case "ban":
                        for ban in objects:
                            channel = guild.get_channel(ban.channel_snowflake)
                            if channel:
                                try:
                                    await channel.set_permissions(
                                        member, overwrite=None
                                    )
                                except discord.Forbidden:
                                    self.__bot.logger.warning(
                                        f"Unable to unban {member.name} ({member.id}) in {channel.name} ({channel.id})."
                                    )
                    case "tmute":
                        for text_mute in objects:
                            channel = guild.get_channel(text_mute.channel_snowflake)
                            if channel:
                                try:
                                    await channel.set_permissions(
                                        member, send_messages=True
                                    )
                                except discord.Forbidden:
                                    self.__bot.logger.warning(
                                        f"Unable to untmute {member.name} ({member.id}) in {channel.name} ({channel.id})."
                                    )
                    case "vmute":
                        for voice_mute in objects:
                            channel = guild.get_channel(voice_mute.channel_snowflake)
                            if channel and member.voice and member.voice.mute:
                                try:
                                    await member.edit(mute=False)
                                except discord.Forbidden:
                                    self.__bot.logger.warning(
                                        f"Unable to unmute {member.name} ({member.id}) in {channel.name} ({channel.id})."
                                    )
                self.__database_factory.delete(**kwargs)

    def __add__(self, pair):
        guild_snowflake, member_snowflake = pair
        self.add_invincible_member(guild_snowflake, member_snowflake)
        return self

    def __sub__(self, pair):
        guild_snowflake, member_snowflake = pair
        self.remove_invincible_member(guild_snowflake, member_snowflake)
        return self

    def add_invincible_member(self, guild_snowflake: int, member_snowflake: int):
        self._invincible_members[(guild_snowflake, member_snowflake)] = True

    @property
    def invincible_members(self):
        return self._invincible_members

    def remove_invincible_member(self, guild_snowflake: int, member_snowflake: int):
        self._invincible_members.pop((guild_snowflake, member_snowflake), None)

    def toggle_enabled(self):
        self.__state = not self.__state
        return self.__state

    def build_clean_dictionary(self):
        dictionary = {}
        for (
            guild_snowflake,
            member_snowflake,
        ), enabled in self._invincible_members.items():
            if not enabled:
                continue
            dictionary.setdefault(guild_snowflake, {"members": {}})
            dictionary[guild_snowflake]["members"][member_snowflake] = True
            skipped_channels = self.__dictionary_service.generate_skipped_members(
                dictionary
            )
            skipped_guilds = self.__dictionary_service.generate_skipped_guilds(
                dictionary
            )
            cleaned_dictionary = self.__dictionary_service.clean_dictionary(
                dictionary=dictionary,
                skipped_channels=skipped_channels,
                skipped_guilds=skipped_guilds,
            )
            # if is_at_home:
            #     if skipped_channels:
            #         pages.extend(
            #             self.__dictionary_service.generate_skipped_dict_pages(
            #                 skipped=skipped_channels,
            #                 title="Skipped Channels in Server",
            #             )
            #         )
            #     if skipped_guilds:
            #         pages.extend(
            #             self.__dictionary_service.generate_skipped_set_pages(
            #                 skipped=skipped_guilds,
            #                 title="Skipped Servers",
            #             )
            #         )
            return cleaned_dictionary

    async def build_pages(self, is_at_home, default_kwargs):
        lines, pages = [], []
        guild_snowflake = default_kwargs.get("guild_snowflake", None)
        guild = self.__bot.get_guild(guild_snowflake)
        title = f"{self.__emoji.get_random_emoji()} {self.__bot.user.display_name} Invincible Members"
        dictionary = self.build_clean_dictionary()
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake in guild_data.get("members", []):
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                lines.append(f"Member: {member.mention}")
                field_count += 1
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
                    )
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)
        return pages
