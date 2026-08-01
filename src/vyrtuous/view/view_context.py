"""!/bin/python3
view_context.py The purpose of this program is to provide context for Views.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord

from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)

INFRACTION_MODELS = [
    ban_service.MODEL,
    flag_service.MODEL,
    text_mute_service.MODEL,
    voice_mute_service.MODEL,
]


class ViewContext:
    def __init__(
        self,
        interaction: discord.Interaction,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
    ):
        self._category: str
        self.channel_snowflake = channel_snowflake
        self.guild_snowflake = guild_snowflake
        self.interaction = interaction
        self.member_snowflake = member_snowflake

    @property
    def category(self):
        return self._category

    @category.setter
    def category(self, new_cat: str):
        if new_cat in [model.identifier for model in INFRACTION_MODELS]:
            self._category = new_cat
        else:
            raise ValueError("Category is not one of ban, flag, tmue or vmute.")
