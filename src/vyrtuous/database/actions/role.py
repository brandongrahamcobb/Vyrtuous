"""flag.py The purpose of this program is to inherit from Action to provide the flag moderation.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.actions.action import Action
from vyrtuous.service.member_snowflake import get_author
from vyrtuous.utils.emojis import get_random_emoji


class Role(Action):

    ACT = "role"
    PLURAL = "Roles"
    SINGULAR = "Role"
    UNDO = "unrole"
    REQUIRED_INSTANTIATION_ARGS = [
        "guild_snowflake",
        "member_snowflake",
        "role_snowflake",
    ]
    OPTIONAL_ARGS = ["channel_snowflake", "created_at", "updated_at"]
    TABLE_NAME = "dummy"

    def __init__(
        self,
        guild_snowflake: int,
        member_snowflake: int,
        role_snowflake: int,
        channel_snowflake: Optional[int] = None,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.role_snowflake = role_snowflake
        self.updated_at = updated_at

    @classmethod
    async def act_embed(cls, action_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(action_information["action_channel_snowflake"])
        author = get_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        role = source.guild.get_member(action_information["action_role_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name} has been granted a role",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, action_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(action_information["action_channel_snowflake"])
        author = get_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        role = source.guild.get_member(action_information["action_role_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name}'s role has been revoked",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
