"""developer.py The purpose of this program is to inherit from the user class as a developer.

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

from datetime import datetime, timezone
from typing import Optional

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.logger import logger
from vyrtuous.utils.emojis import get_random_emoji

class Bug(DatabaseFactory):

    ACT = "bug"
    CATEGORY = "bug"
    PLURAL = "Bugs"
    SCOPES = ["channels"]
    SINGULAR = "Bug"
    UNDO = None
    
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "id",
        "member_snowflakes",
        "message_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "notes",
        "resolved",
        "updated_at",
    ]
    TABLE_NAME = "bug_tracking"

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        id: str,
        member_snowflakes: list[int],
        message_snowflake: int,
        created_at: Optional[datetime] = None,
        notes: Optional[str] = None,
        resolved: Optional[bool] = False,
        updated_at: Optional[datetime] = None,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.created_at: datetime = created_at or datetime.now(timezone.utc)
        self.member_snowflakes = member_snowflakes
        self.id = id
        self.guild_snowflake = guild_snowflake
        self.message_snowflake = message_snowflake
        self.notes = notes
        self.resolved = resolved
        self.updated_at: datetime = updated_at or datetime.now(timezone.utc)

    async def create_embed(self, action, member_snowflake, source):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(self.channel_snowflake)
        guild = bot.get_guild(self.guild_snowflake)
        try:
            msg = await channel.fetch_message(self.message_snowflake)
        except discord.NotFound:
            logger.warning(f"Message reference not found ({self.message_snowflake}).")
        author = resolve_author(source=source)
        member = bot.get_user(member_snowflake)
        user_mentions = []
        for member_snowflake in self.member_snowflakes:
            user_mentions.append(bot.get_user(member_snowflake).mention)
        embed = discord.Embed(
            title=f"{get_random_emoji()} {member.display_name} has been {action}",
            description=(
                f"**Guild:** {guild.name}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Message:** {msg.jump_url}\n"
                f"**By:** {author.mention}\n"
                f"**Assigned devs:** {', '.join(user_mentions)}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed