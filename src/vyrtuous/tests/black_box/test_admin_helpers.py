
''' test_admin_helpers.py The purpose of this program is to provide the shared test variables for tests using Discord objects.
    Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
'''
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot

async def admin_cleanup(guild_id: Optional[int], privileged_author_id: Optional[int]):
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        await conn.execute('''
            UPDATE users SET administrator_guild_ids=$2 WHERE discord_snowflake=$1
        ''', int(privileged_author_id), [int(guild_id)])

async def admin_initiation(guild_id: Optional[int], privileged_author_id: Optional[int]):
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        await conn.execute('''
            INSERT INTO users (discord_snowflake, developer_guild_ids, updated_at, created_at)
            VALUES ($1, $2, NOW(), NOW())
            ON CONFLICT (discord_snowflake) 
            DO UPDATE SET developer_guild_ids = $2, updated_at = NOW()
        ''', int(privileged_author_id), [int(guild_id)])