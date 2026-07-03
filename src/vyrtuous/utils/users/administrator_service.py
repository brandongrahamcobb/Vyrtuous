"""!/bin/python3
administrator_service.py The purpose of this program is to extend Service to service the administrator and administrator role classes.

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
from vyrtuous.db.administrator import Administrator, NotAdministrator
from vyrtuous.db.database_factory import DatabaseFactory

MODEL = Administrator


async def is_administrator(guild_snowflake: int, member_snowflake: int) -> bool:
    database_factory = DatabaseFactory(MODEL)
    administrator = await database_factory.select(
        guild_snowflake=int(guild_snowflake),
        member_snowflake=int(member_snowflake),
        singular=True,
    )
    if not administrator:
        raise NotAdministrator
    return True


async def is_administrator_wrapper(
    context,
):
    return await is_administrator(
        guild_snowflake=int(context.guild_snowflake),
        member_snowflake=int(context.member_snowflake),
    )


async def is_administrator_at_all(
    member_snowflake: int,
):
    database_factory = DatabaseFactory(MODEL)
    administrator = await database_factory.select(
        member_snowflake=int(member_snowflake), singular=True
    )
    if not administrator:
        raise NotAdministrator
    return True


async def administrators_by_role(role_snowflake: int):
    database_factory = DatabaseFactory(MODEL)
    administrators = await database_factory.select(
        role_snowflakes=int(role_snowflake),
        inside_fields=["role_snowflakes"],
        singular=False,
    )
    if administrators:
        return administrators
    return []


async def administrator_existing(
    guild_snowflake: int, member_snowflake: int, role_snowflake: int
):
    database_factory = DatabaseFactory(MODEL)
    administrator = await database_factory.select(
        guild_snowflake=int(guild_snowflake),
        member_snowflake=int(member_snowflake),
        role_snowflakes=[int(role_snowflake)],
        inside_fields=["role_snowflakes"],
        singular=True,
    )
    if administrator:
        return administrator
    return None


async def added_role(guild_snowflake: int, member_snowflake: int, role_snowflake: int):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    administrator_role_snowflakes = []
    administrator = await administrator_existing(
        guild_snowflake=guild_snowflake,
        member_snowflake=member_snowflake,
        role_snowflake=role_snowflake,
    )
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        return
    if not administrator:
        administrator = Administrator(
            guild_snowflake=int(guild_snowflake),
            member_snowflake=int(member_snowflake),
            role_snowflakes=[int(role_snowflake)],
        )
        await database_factory.create(administrator)
        return
    administrator_role_snowflakes = administrator.role_snowflakes
    administrator_role_snowflakes.append(role_snowflake)
    where_kwargs = {
        "guild_snowflake": int(guild_snowflake),
        "member_snowflake": member_snowflake,
    }
    set_kwargs = {"role_snowflakes": administrator_role_snowflakes}
    await database_factory.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
    bot.logger.info(
        f"Granted administrator to member ({member_snowflake}) in guild ({guild_snowflake})."
    )


async def removed_role(
    guild_snowflake: int, member_snowflake: int, role_snowflake: int
):
    bot = DiscordBot.get_instance()
    database_factory = DatabaseFactory(MODEL)
    administrator_role_snowflakes = []
    administrator = await database_factory.select(
        guild_snowflake=int(guild_snowflake),
        member_snowflake=int(member_snowflake),
        role_snowflakes=role_snowflake,
        singular=True,
        inside_fields=["role_snowflakes"],
    )
    if not administrator:
        return
    administrator_role_snowflakes = administrator.role_snowflakes
    administrator_role_snowflakes.remove(role_snowflake)
    if administrator_role_snowflakes == []:
        await database_factory.delete(
            guild_snowflake=guild_snowflake, member_snowflake=member_snowflake
        )
    else:
        where_kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "member_snowflake": member_snowflake,
        }
        set_kwargs = {"role_snowflakes": administrator_role_snowflakes}
        await database_factory.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
    bot.logger.info(
        f"Revoked administrator from member ({member_snowflake}) in guild ({guild_snowflake})."
    )
