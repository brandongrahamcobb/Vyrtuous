"""!/bin/python3stage"
clean_automute_service.py The purpose of this program is to clean expired automutes.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.video_channel import VideoChannel

MODEL = VideoChannel


async def clean_expired_video_only_channels() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    expired_video_only_channels = await database_factory.select(
        expired=True, singular=False
    )
    if expired_video_only_channels:
        for expired_video_only_channel in expired_video_only_channels:
            await database_factory.delete(
                channel_snowflake=expired_video_only_channel.channel_snowflake,
                guild_snowflake=expired_video_only_channel.guild_snowflake,
            )
        bot.logger.info("Cleaned up expired video-only channels.")
    return count
