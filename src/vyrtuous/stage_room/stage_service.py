from copy import copy

"""!/bin/python3
stage_service.py The purpose of this program is to extend Service to service the stage class.

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

import time

import discord

from vyrtuous.stage_room.stage import Stage

# from vyrtuous.voice_mute.voice_mute import VoiceMute


class StageService:
    __CHUNK_SIZE = 7
    MODEL = Stage

    def __init__(
        self,
        *,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        duration=None,
        emoji=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__duration = duration
        self.__emoji = emoji

    async def send_stage_ask_to_speak_message(
        self, join_log: dict[int, discord.Member], member: discord.Member, stage: Stage
    ):
        now = time.time()
        join_log[member.id] = [t for t in join_log[member.id] if now - t < 300]
        if len(join_log[member.id]) < 1:
            join_log[member.id].append(now)
            embed = discord.Embed(
                title=f"{self.__emoji.get_random_emoji()} {stage.channel_snowflake} — Stage Mode",
                description=f"Ends <t:{int(stage.expires_in.timestamp())}:R>",
                color=discord.Color.green(),
            )
            embed.add_field(name="\u200b", value="**Ask to speak!**", inline=False)
            await self.__bot.get_channel(stage.channel_snowflake).send(embed=embed)

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        pages = []
        dictionary = {}
        stages = await self.__database_factory.select(singular=False, **where_kwargs)
        for stage in stages:
            dictionary.setdefault(stage.guild_snowflake, {"channels": {}})
            dictionary[stage.guild_snowflake]["channels"].setdefault(
                stage.channel_snowflake, {}
            )
            dictionary[stage.guild_snowflake]["channels"][
                stage.channel_snowflake
            ].setdefault("stages", {})
            dictionary[stage.guild_snowflake]["channels"][stage.channel_snowflake][
                "stages"
            ].update({"expires_in": self.__duration.from_expires_in(stage.expires_in)})
        skipped_channels = self.__dictionary_service.generate_skipped_channels(
            dictionary
        )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        # if is_at_home:
        #     if skipped_guilds:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_set_pages(
        #                 skipped=skipped_guilds,
        #                 title="Skipped Servers",
        #             )
        #         )
        #     if skipped_channels:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_dict_pages(
        #                 skipped=skipped_channels,
        #                 title="Skipped Channels in Server",
        #             )
        #         )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Stages"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )
        embed = discord.Embed(
            title=title, description="Default view", color=discord.Color.blue()
        )
        if not dictionary:
            pages = [embed]
        stage_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, stage_dictionary in guild_data.get(
                "channels"
            ).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                lines.append(
                    f"**Expires in:** {stage_dictionary.get('stages', {}).get('expires_in', None)}"
                )
                stage_n += 1
                field_count += 1
                if field_count == self.__CHUNK_SIZE:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = flush_page(embed, pages, title, guild.name)
                    lines = []
                    field_count = 0
                if lines:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({stage_n})**"
        return pages

    async def toggle_stage(self, channel_dict, default_kwargs, duration):
        updated_kwargs = default_kwargs.copy()
        guild = self.__bot.get_guild(updated_kwargs.get("guild_snowflake", None))
        failed, pages, skipped, succeeded = [], [], [], []
        updated_kwargs.update(channel_dict.get("columns", None))
        stage_kwargs = updated_kwargs.copy()
        del stage_kwargs["member_snowflake"]
        stage = await self.__database_factory.select(**stage_kwargs, singular=True)
        if stage:
            title = f"{self.__emoji.get_random_emoji()} Stage Ended in {channel_dict.get('mention', None)}"
            await self.__database_factory.delete(**updated_kwargs)
            # for member in channel_dict.get("object", None).members:
            # await VoiceMute.delete(
            #     **updated_kwargs,
            #     member_snowflake=member.id,
            #     target="room",
            # )
            # voice_mute = await VoiceMute.select(
            #     **updated_kwargs,
            #     member_snowflake=member.id,
            #     target="user",
            #     singular=True,
            # )
            # if not voice_mute and member.voice and member.voice.mute:
            #     try:
            #         await member.edit(
            #             mute=False,
            #             reason="Stage ended — no user-specific mute found",
            #         )
            #         succeeded.append(member)
            #     except discord.Forbidden as e:
            #         logger.warning(
            #             f"Unable to undo voice-mute "
            #             f"for member {member.display_name} ({member.id}) in "
            #             f"channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) in "
            #             f"guild {guild.name} ({guild.id}). "
            #             f"{str(e).capitalize()}"
            #         )
            #         failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict.get('mention', None)}",
                f"**Unmuted:** {len(succeeded)} users",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=title,
                color=discord.Color.blurple(),
            )
            pages.append(embed)
        else:
            stage = self.MODEL(
                **stage_kwargs,
                expires_in=duration.expires_in,
            )
            await self.__database_factory.create(stage)
            # for member in channel_dict.get("object", None).members:
            # if await PermissionService.check(
            #     **updated_kwargs,
            #     lowest_role="Coordinator",
            # ) or member.id == updated_kwargs.get("member_snowflake", None):
            #     skipped.append(member)
            #     continue
            # voice_mute = await VoiceMute(
            #     **updated_kwargs,
            #     expires_in=duration.expires_in,
            #     member_snowflake=member.id,
            #     target="room",
            #     reason="Stage mute",
            # )
            # await voice_mute.create()
            # try:
            #     if member.voice and member.voice.channel.id == channel_dict.get(
            #         "id", None
            #     ):
            #         await member.edit(mute=True)
            #     succeeded.append(member)
            # except Exception as e:
            #     logger.warning(
            #         f"Unable to voice-mute "
            #         f"member {member.display_name} ({member.id}) "
            #         f"in channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) "
            #         f"in guild {guild.name} ({guild.id}). "
            #         f"{str(e).capitalize()}"
            #     )
            #     failed.append(member)
            description_lines = [
                f"**Channel:** {channel_dict.get('mention', None)}",
                f"**Expires:** {duration}",
                f"**Muted:** {len(succeeded)} users",
                f"**Skipped:** {len(skipped)}",
            ]
            if failed:
                description_lines.append(f"**Failed:** {len(failed)}")
            embed = discord.Embed(
                description="\n".join(description_lines),
                title=f"{self.__emoji.get_random_emoji()} Stage Created in {channel_dict.get('name', None)}",
                color=discord.Color.blurple(),
            )
            pages.append(embed)
        return pages

    async def toggle_stage_mute(self, channel_dict, default_kwargs, member_dict):
        # await PermissionService.has_equal_or_lower_role(
        #     **default_kwargs,
        #     target_member_snowflake=member_dict.get("id", None),
        # )
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        stage = await self.__database_factory.select(singular=True, **updated_kwargs)
        if stage:
            await member_dict.get("object", None).edit(
                mute=not member_dict.get("object", None).voice.mute
            )
            return f"Successfully toggled the mute for {member_dict.get('mention', None)} in {channel_dict.get('mention', None)}."
