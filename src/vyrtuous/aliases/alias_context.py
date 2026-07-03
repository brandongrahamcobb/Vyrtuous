"""!/bin/python3
alias_context.py The purpose of this program is to classify Alias context.

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

from typing import Union

from vyrtuous.aliases import alias_service
from vyrtuous.aliases.ban_alias import BanAlias
from vyrtuous.aliases.flag_alias import FlagAlias
from vyrtuous.aliases.role_alias import RoleAlias
from vyrtuous.aliases.text_mute_alias import TextMuteAlias
from vyrtuous.aliases.voice_mute_alias import VoiceMuteAlias
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.alias import Alias, NotAlias
from vyrtuous.db.database_factory import DatabaseFactory

MODEL = Alias


class AliasContext:

    def __init__(self, content: str, guild_snowflake: int):
        self.__alias_name: str
        self.__content = content
        self.__args = []
        self.__kwargs: dict[str, tuple] = {}
        self.category: str
        self.channel_snowflake: int
        self.duration_value: str
        self.guild_snowflake: int = guild_snowflake
        self.member_snowflake: int
        self.role_snowflake: int
        self.reason = "No reason provided"

    async def setup(self):
        self._message_to_args()
        self._alias_name_from_args()
        alias = await self._populate_alias()
        if not alias:
            return False
        self._fill_map(alias)
        await self._convert_args_to_values()
        return True

    def _message_to_args(self) -> None:
        bot = DiscordBot.get_instance()
        self.__args = (
            self.__content[len(bot.config["discord_command_prefix"]) :].strip().split()
        )

    def _fill_map(
        self,
        alias: Union[BanAlias, FlagAlias, RoleAlias, TextMuteAlias, VoiceMuteAlias],
    ) -> None:
        map = alias.ARGS_MAP
        sorted_args = sorted(map.items(), key=lambda x: x[1])
        for i, (key, pos) in enumerate(sorted_args):
            if i == len(sorted_args) - 1:
                value = (
                    " ".join(str(a) for a in self.__args[pos - 1 :])
                    if len(self.__args) >= pos
                    else ""
                )
            else:
                value = str(self.__args[pos - 1]) if len(self.__args) >= pos else ""
            self.__kwargs[key] = (pos, value)

    def _alias_name_from_args(self):
        if not self.__args:
            return
        self.__alias_name = self.__args[0]

    async def _populate_alias(
        self,
    ) -> Union[BanAlias, FlagAlias, RoleAlias, TextMuteAlias, VoiceMuteAlias]:
        database_factory = DatabaseFactory(MODEL)
        alias_entry = await database_factory.select(
            alias_name=self.__alias_name,
            guild_snowflake=self.guild_snowflake,
            singular=True,
        )
        if not alias_entry:
            raise NotAlias
        self.category = alias_entry.category
        self.channel_snowflake = alias_entry.channel_snowflake
        self.role_snowflake = alias_entry.role_snowflake
        alias = alias_service.alias_category_to_alias(category=alias_entry.category)
        return alias

    async def _convert_args_to_values(self):
        for field, tuple in self.__kwargs.items():
            value = tuple[1]
            if field == "duration":
                if not value:
                    self.duration_value = "8h"
                else:
                    self.duration_value = value
            elif field == "member":
                self.member_snowflake = int(value.replace("<@", "").replace(">", ""))
            elif field == "reason":
                if value is None:
                    self.reason = "No reason provided."
                else:
                    self.reason = value
