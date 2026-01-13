"""action.py The purpose of this program is to be a child of DatabaseFactory and the parent to all moderation actions.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
from typing import Optional

from vyrtuous.database.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot

class Action(DatabaseFactory):

    ACT = "alias"
    PLURAL = "Aliases"
    SINGULAR = "Alias"
    UNDO = "alias"
    REQUIRED_INSTANTIATION_ARGS = [
        "alias_name",
        "alias_type"
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

    def __init__(
        self,
        alias_name: Optional[str],
        alias_type: Optional[str],
        channel_snowflake: Optional[int],
        guild_snowflake: Optional[int],
        role_snowflake: Optional[int] = None,
    ):
        self.alias_cog = self.bot.get_cog("Aliases")
        self.alias_type = alias_type
        self.alias_name = alias_name
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>"
        self.guild_snowflake = guild_snowflake
        self.handlers = {
            "ban": self.alias_cog.handle_ban_alias,
            "vegan": self.alias_cog.handle_vegan_alias,
            "carnist": self.alias_cog.handle_carnist_alias,
            "unban": self.alias_cog.handle_unban_alias,
            "flag": self.alias_cog.handle_flag_alias,
            "unflag": self.alias_cog.handle_unflag_alias,
            "voice_mute": self.alias_cog.handle_voice_mute_alias,
            "unvoice_mute": self.alias_cog.handle_unmute_alias,
            "text_mute": self.alias_cog.handle_text_mute_alias,
            "untext_mute": self.alias_cog.handle_untextmute_alias,
            "role": self.alias_cog.handle_role_alias,
            "unrole": self.alias_cog.handle_unrole_alias,
        }
        self.handler = self.handlers[alias_type]
        self.role_snowflake = role_snowflake
        self.role_mention = f"<@&{role_snowflake}>"

    @property
    def alias_type(self):
        return self._alias_type

    @alias_type.setter
    def alias_type(self, alias_type: Optional[str]):
        if alias_type not in (
            "vegan",
            "carnist",
            "voice_mute",
            "unvoice_mute",
            "ban",
            "unban",
            "flag",
            "unflag",
            "text_mute",
            "untext_mute",
            "role",
            "unrole",
        ):
            raise ValueError("Invalid alias_type.")
        self._alias_type = alias_type
        
    @classmethod
    def format_aliases(cls, aliases) -> list[str]:
        lines = []
        if not aliases:
            return []
        grouped = defaultdict(list)
        for alias in aliases:
            match alias.alias_type:
                case "ban" | "unban":
                    formatted_type = "Ban"
                case "vegan" | "carnist":
                    formatted_type = "Veganism"
                case "role" | "unrole":
                    formatted_type = "Role"
                case "flag" | "unflag":
                    formatted_type = "Flag"
                case "text_mute" | "untext_mute":
                    formatted_type = "Text Mute"
                case "voice_mute" | "unvoice_mute":
                    formatted_type = "Voice Mute"
            grouped[(alias.channel_snowflake, formatted_type)].append(alias)
        for (channel_snowflake, formatted_type), channel_aliases in grouped.items():
            lines.append(f"**{formatted_type}**")
            for alias in sorted(channel_aliases, key=lambda a: a.alias_name.lower()):
                if alias.role_snowflake:
                    lines.append(f"`{alias.alias_name}` → <@&{alias.role_snowflake}>")
                else:
                    lines.append(f"`{alias.alias_name}`")
        return lines
    