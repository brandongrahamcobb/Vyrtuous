from datetime import datetime
from typing import Union

import discord

from vyrtuous.ban.ban import Ban
from vyrtuous.base.context import Context
from vyrtuous.flag.flag import Flag
from vyrtuous.text_mute.text_mute import TextMute
from vyrtuous.utils.permission_service import PermissionService
from vyrtuous.utils.view_utilities import limit_available_to_top_25_by_member_count
from vyrtuous.voice_mute.voice_mute import VoiceMute


class ViewContext(Context):
    def __init__(self, interaction):
        super().__init__(interaction=interaction)
        self.available_channels: list[discord.VoiceChannel | None] = []
        self.available_guilds: list[discord.Guild | None] = []
        self._expires_in: datetime | None = None
        self.infraction = None
        self.interaction = interaction
        self._record: Union[Ban, Flag, TextMute, VoiceMute] | None = None
        self._target_member_snowflake: int | None = None
        self._target_channel_snowflake: int | None = None
        self.service = None

    async def setup(self, target_member_snowflake):
        self.target_member_snowflake = target_member_snowflake
        available_channels, available_guilds = await PermissionService.can_list(
            source=self.interaction
        )
        self.available_channels = limit_available_to_top_25_by_member_count(
            available=available_channels
        )
        self.available_guilds = limit_available_to_top_25_by_member_count(
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
    def record(self):
        return self._record

    @record.setter
    def record(self, nr):
        if isinstance(nr, (Ban, Flag, TextMute, VoiceMute)):
            self._record = nr
        else:
            raise ValueError("Record is not one of Ban, Flag, TextMute or VoiceMute.")
