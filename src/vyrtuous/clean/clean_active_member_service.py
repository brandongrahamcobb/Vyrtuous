"""!/bin/python3
active_member.py The purpose of this program is to service active members.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from datetime import datetime, timedelta, timezone

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.active_member import ActiveMember
from vyrtuous.db.database_factory import DatabaseFactory

MODEL = ActiveMember


async def clean_inactive_members() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    saved_members = await database_factory.select(singular=False)
    for member in saved_members:
        if datetime.now(timezone.utc) - member.last_active > timedelta(days=7):
            del bot.registry.get(MemberState).active[member.member_snowflake]
            count += 1
            bot.logger.info(
                f"Deleted inactive member {member.display_name} from the active members database table."
            )
    return count
