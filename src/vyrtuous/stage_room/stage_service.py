"""!/bin/python3stage"
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
from collections import defaultdict
from copy import copy
from dataclasses import dataclass, field
from typing import Any, Dict, List

import discord

from vyrtuous.stage_room.stage import Stage


@dataclass
class StageDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[str, Any]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


class StageService:
    __CHUNK_SIZE = 7
    MODEL = Stage

    def __init__(
        self,
        *,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        duration_service=None,
        emoji=None,
        moderator_service=None,
        voice_mute_service=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__duration_service = duration_service
        self.__emoji = emoji
        self.__moderator_service = moderator_service
        self.__voice_mute_service = voice_mute_service
        self.__join_log = defaultdict(list)

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

    async def build_dictionary(self, where_kwargs):
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
            ].update(
                {
                    "expires_in": self.__duration_service.from_expires_in(
                        stage.expires_in
                    )
                }
            )
        return dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Stages"

        where_kwargs = object_dict.get("columns", None)

        dictionary = await self.build_dictionary(where_kwargs=where_kwargs)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=StageDictionary, dictionary=dictionary
        )
        if is_at_home:
            pages.extend(processed_dictionary.skipped_channels)
            pages.extend(processed_dictionary.skipped_guilds)

        stage_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
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
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
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

    async def toggle_stage(self, channel_dict, context, duration):
        failed, pages, skipped, succeeded = [], [], [], []
        stage = await self.__database_factory.select(
            **channel_dict.get("columns", None), singular=True
        )
        if stage:
            title = f"{self.__emoji.get_random_emoji()} Stage Ended in {channel_dict.get('mention', None)}"
            await self.__database_factory.delete(**channel_dict.get("columns", None))
            failed, succeeded = await self.__voice_mute_service.off_stage(
                channel_dict=channel_dict
            )
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
                **channel_dict.get("columns", None),
                expires_in=self.__duration_service.to_expires_in(duration),
            )
            await self.__database_factory.create(stage)
            failed, skipped, succeeded = await self.__voice_mute_service.on_stage(
                channel_dict=channel_dict,
                context=context,
                duration=duration,
            )
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

    async def toggle_stage_mute(self, channel_dict, context, member_dict):
        await self.__moderator_service.has_equal_or_lower_role(
            **context.to_dict(),
            target_member_snowflake=member_dict.get("id", None),
        )
        stage = await self.__database_factory.select(
            singular=True, **channel_dict.get("columns", None)
        )
        if stage:
            await member_dict.get("object", None).edit(
                mute=not member_dict.get("object", None).voice.mute
            )
            return f"Successfully toggled the mute for {member_dict.get('mention', None)} in {channel_dict.get('mention', None)}."

    async def clean_expired(self):
        expired_stages = await self.__database_factory.select(expired=True)
        if expired_stages:
            for expired_stage in expired_stages:
                channel_snowflake = int(expired_stage.channel_snowflake)
                guild_snowflake = int(expired_stage.guild_snowflake)
                guild = self.__bot.get_guild(guild_snowflake)
                if guild is None:
                    await self.__database_factory.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                    )
                    self.__bot.logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired stage."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await self.__database_factory.delete(
                        channel_snowflake=channel_snowflake,
                        guild_snowflake=guild_snowflake,
                    )
                    self.__bot.logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                    )
                    continue
                self.__voice_mute_service.clean_stage_expired(
                    channel=channel, guild=guild
                )
            self.__bot.logger.info("Cleaned up expired stages.")

    async def migrate(self, kwargs):
        self.__database_factory.update(**kwargs)

    async def enforce(self, after, member):
        should_be_muted = False
        expires_in = None
        stage = await self.__database_factory.select(channel_snowlfake=after.channel)
        if stage:
            await self.send_stage_ask_to_speak_message(
                join_log=self.__join_log, member=member, stage=stage
            )
            kwargs = {
                "channel_snowflake": after.channel.id,
                "guild_snowflake": after.channel.guild.id,
                "member_snowflake": member.id,
            }
            highest_role = await self.__moderator_service.resolve_highest_role(
                kwargs=kwargs
            )
            if highest_role == "Everyone":
                should_be_muted = True
                expires_in = stage.expires_in
        return should_be_muted, expires_in

    async def is_active_stage_room(self, channel):
        stage_room = await self.__database_factory.select(channel_snowflake=channel.id)
        if stage_room:
            return True
        return False
