"""vegan.py The purpose of this program is to inherit from the DatabaseFactory to provide the vegan role.

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
from typing import Optional

import discord

from vyrtuous.database.database_factory import DatabaseFactory
from vyrtuous.service.member_snowflake import get_author


class Vegan(DatabaseFactory):

    ACT = "vegan"
    PLURAL = "Vegans"
    SINGULAR = "Vegan"
    UNDO = "carnist"
    REQUIRED_INSTANTIATION_ARGS = ["guild_snowflake", "member_snowflake"]
    OPTIONAL_ARGS = ["created_at", "updated_at"]
    TABLE_NAME = "vegans"

    def __init__(
        self,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.updated_at = updated_at

    @classmethod
    async def act_embed(cls, action_information, source, **kwargs):
        author = get_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        embed = discord.Embed(
            title=f"\U0001f525\U0001f525 {member.display_name} "
            f"is going Vegan!!!\U0001f525\U0001f525",
            description=(f"**By:** {author.mention}\n" f"**User:** {member.mention}\n"),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, action_information, source, **kwargs):
        author = get_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        embed = discord.Embed(
            title=f"\U0001f44e\U0001f44e "
            f"{member.display_name} is a Carnist \U0001f44e\U0001f44e",
            description=(f"**By:** {author.mention}\n" f"**User:** {member.mention}\n"),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
