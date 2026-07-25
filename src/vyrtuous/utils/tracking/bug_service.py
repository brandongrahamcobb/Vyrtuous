"""!/bin/python3
bug_service.py The purpose of this program is to service bugs.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.bug import Bug
from vyrtuous.db.database_factory import DatabaseFactory

MODEL = Bug


async def create_bug(message, reference) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    try:
        bug = MODEL(
            channel_snowflake=message.channel.id,
            member_snowflakes=[],
            guild_snowflake=message.guild.id,
            id=reference,
            message_snowflake=message.id,
        )
        await database_factory.create(bug)
    except discord.Forbidden as e:
        bot.logger.info(str(e).capitalize())
