"""!/bin/python3

view_context.py The purpose of this program is to provide context for Views.

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

from datetime import datetime
from typing import Union

import discord

from vyrtuous.base.context import Context


class ViewContext(Context):
    def __init__(
        self,
        *,
        ban_service=None,
        flag_service=None,
        interaction: discord.Interaction = None,
        moderator_service=None,
        text_mute_service=None,
        voice_mute_service=None,
    ):
        super().__init__(interaction=interaction)
        self.available_channels: list[discord.VoiceChannel | None] = []
        self.available_guilds: list[discord.Guild | None] = []
        self._expires_in: datetime | None = None
        self.interaction = interaction
        self._target_member_snowflake: int | None = None
        self._target_channel_snowflake: int | None = None
        self.service = None
        self.__ban_service = ban_service
        self.__voice_mute_service = voice_mute_service
        self.__flag_service = flag_service
        self.__text_mute_service = text_mute_service
        self._infraction = Union[
            self.__ban_service.MODEL,
            self.__flag_service.MODEL,
            self.__text_mute_service.MODEL,
            self.__voice_mute_service.MODEL,
        ]
        self.__moderator_service = moderator_service

    def limit_available_to_top_25_by_member_count(self, available):
        all_key = "all"
        items = []
        if all_key in available:
            items.append(all_key)
        objects = []
        for k, v in available.items():
            if k == all_key:
                continue
            if isinstance(v, (list, tuple, set)):
                objects.extend(v)
            else:
                objects.append(v)
        objects.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
        top_25 = objects[:25]
        return items + top_25

    async def setup(self, target_member_snowflake):
        self.target_member_snowflake = target_member_snowflake
        available_channels, available_guilds = await self.__moderator_service.can_list(
            source=self.interaction
        )
        self.available_channels = self.limit_available_to_top_25_by_member_count(
            available=available_channels
        )
        self.available_guilds = self.limit_available_to_top_25_by_member_count(
            available=available_guilds
        )

    @property
    def target_member_snowflake(self):
        return self._target_member_snowflake

    @target_member_snowflake.setter
    def target_member_snowflake(self, ntms):
        if isinstance(ntms, int):
            self._target_member_snowflake = ntms
        else:
            raise ValueError("Target member snowflake is not an integer.")

    @property
    def target_channel_snowflake(self):
        return self._target_channel_snowflake

    @target_channel_snowflake.setter
    def target_channel_snowflake(self, ntcs):
        if isinstance(ntcs, int):
            self._target_channel_snowflake = ntcs
        else:
            raise ValueError("Target channel snowflake is not an integer.")

    @property
    def expires_in(self):
        return self._expires_in

    @expires_in.setter
    def expires_in(self, neo):
        if isinstance(neo, datetime):
            self._expires_in = neo
        else:
            raise ValueError("Expires in is not an datetime.")

    @property
    def infraction(self):
        return self._infraction

    @infraction.setter
    def infraction(self, nr):
        if isinstance(
            nr,
            (
                self.__ban_service.MODEL,
                self.__flag_service.MODEL,
                self.__text_mute_service.MODEL,
                self.__voice_mute_service.MODEL,
            ),
        ):
            self._record = nr
        else:
            raise ValueError(
                "Infraction is not one of Ban, Flag, TextMute or VoiceMute."
            )
