"""!/bin/python3
active_member.py The purpose of this program is to service active members.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.active_member import ActiveMember
from vyrtuous.db.database_factory import DatabaseFactory

MODEL = ActiveMember


async def is_active(member_snowflake: int) -> bool:
    bot: DiscordBot = DiscordBot.get_instance()
    if member_snowflake in bot.registry.get(MemberState).active.keys():
        return True
    return False


async def populate() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    active_members = await database_factory.select(singular=False)
    for member in active_members:
        bot.registry.get(MemberState).active[
            member.member_snowflake
        ] = member.display_name
    bot.logger.info("Populated in-memory active members.")


async def save_and_update_active_members() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    saved_members = await database_factory.select(singular=False)
    member_snowflakes = [
        active_member.member_snowflake for active_member in saved_members
    ]
    for member_snowflake, data in bot.registry.get(MemberState).active.items():
        if member_snowflake not in member_snowflakes:
            active_member = ActiveMember(
                display_name=data[0],
                last_active=data[1],
                member_snowflake=int(member_snowflake),
            )
            await database_factory.create(active_member)
            bot.logger.info(f"Saved {data[0]} to active members database table.")
        else:
            where_kwargs = {
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"last_active": data[1]}
            await database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )
            bot.logger.info(
                f"Updated {data[0]} last active timestamp in the active members database table."
            )
