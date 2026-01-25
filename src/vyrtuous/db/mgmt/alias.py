"""action.py The purpose of this program is to be a child of DatabaseFactory and the parent to all moderation actions.

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

from collections import defaultdict
from datetime import datetime, timezone
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.actions.ban import Ban
from vyrtuous.db.actions.flag import Flag
from vyrtuous.db.actions.hide import Hide
from vyrtuous.db.actions.role import Role
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.fields.duration import DurationObject


class Alias(DatabaseFactory):

    ACT = "alias"
    CATEGORY = "alias"
    PLURAL = "Aliases"
    SINGULAR = "Alias"
    SCOPES = ["channels"]
    UNDO = "xalias"

    REQUIRED_INSTANTIATION_ARGS = [
        "alias_name",
        "alias_type",
        "channel_snowflake",
        "guild_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "expired",
        "expires_in",
        "reason",
        "role_snowflake",
        "updated_at",
    ]

    TABLE_NAME = "command_aliases"

    _ALIAS_CLASS_MAP = {
        "ban": Ban,
        "hide": Hide,
        "vmute": VoiceMute,
        "tmute": TextMute,
        "role": Role,
        "flag": Flag,
        "vegan": Vegan,
    }

    def __init__(
        self,
        alias_name: str,
        alias_type: str,
        channel_snowflake: int,
        guild_snowflake: int,
        created_at: Optional[datetime] = None,
        role_snowflake: Optional[int] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.bot = DiscordBot.get_instance()
        self.alias_class = self._ALIAS_CLASS_MAP.get(alias_type, None)
        self.alias_cog = self.bot.get_cog("Aliases")
        self.alias_type = alias_type
        self.alias_name = alias_name
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>"
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.handlers = {
            "ban": self.alias_cog.handle_ban_alias,
            "vegan": self.alias_cog.handle_vegan_alias,
            "carnist": self.alias_cog.handle_carnist_alias,
            "unban": self.alias_cog.handle_unban_alias,
            "flag": self.alias_cog.handle_flag_alias,
            "hide": self.alias_cog.handle_hide_alias,
            "unhide": self.alias_cog.handle_unhide_alias,
            "unflag": self.alias_cog.handle_unflag_alias,
            "vmute": self.alias_cog.handle_voice_mute_alias,
            "unvmute": self.alias_cog.handle_unmute_alias,
            "tmute": self.alias_cog.handle_text_mute_alias,
            "untmute": self.alias_cog.handle_untextmute_alias,
            "role": self.alias_cog.handle_role_alias,
            "unrole": self.alias_cog.handle_unrole_alias,
        }
        self.role_snowflake = role_snowflake
        self.role_mention = f"<@&{role_snowflake}>"
        self.updated_at = updated_at

    @property
    def alias_type(self):
        return self._alias_type

    @alias_type.setter
    def alias_type(self, alias_type: str):
        if alias_type not in (
            "vegan",
            "vmute",
            "ban",
            "hide",
            "flag",
            "tmute",
            "role",
        ):
            raise ValueError("Invalid alias_type.")
        self._alias_type = alias_type

    @classmethod
    async def generate_cap_duration(
        cls,
        channel_snowflake: int,
        guild_snowflake: int,
        moderation_type: str,
    ):
        cap = await Cap.select(
            channel_snowflake=channel_snowflake,
            guild_snowflake=guild_snowflake,
            moderation_type=moderation_type,
            singular=True,
        )
        if not hasattr(cap, "duration"):
            cap_duration = DurationObject("8h").to_seconds()
        else:
            cap_duration = cap.duration_seconds
        return cap_duration