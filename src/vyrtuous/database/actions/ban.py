"""ban.py The purpose of this program is to inherit from Action to provide the ban moderation.

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


class Ban(Action):

    ACT = "ban"
    PLURAL = "Bans"
    SINGULAR = "Ban"
    UNDO = "unban"
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "member_snowflake",
    ]
    OPTIONAL_ARGS = ["created_at", "expires_in", "reason", "updated_at"]
    TABLE_NAME = "active_bans"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: Optional[datetime] = None,
        expired: bool = False,
        expires_in: Optional[datetime] = None,
        reason: Optional[str] = "No reason provided.",
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.expired = expired
        self.expires_in = expires_in
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.reason = reason
        self.updated_at = updated_at

    async def create_embed(self, action_information, source, **kwargs):
        channel = self.bot.get_channel(action_information["action_channel_snowflake"])
        author = get_author(source=source)
        member = self.bot.get_member(action_information["action_member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Banned",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {action_information['action_duration']}\n"
                f"**Reason:** {action_information['action_reason']}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    def get_handler(cls):
        bot = DiscordBot.get_instance()
        alias_cog = bot.get_cog("Aliases")
        return alias_cog.handle_ban_alias
