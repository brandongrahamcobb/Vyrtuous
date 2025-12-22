
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
from vyrtuous.inc.helpers import *

async def coord_cleanup(channel_id: Optional[int], guild_id: Optional[int], privileged_author_id: Optional[int]):
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        await conn.execute('''
            UPDATE users SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $2), updated_at = NOW() WHERE discord_snowflake = $1
        ''', privileged_author_id, channel_id)

async def coord_initiation(channel_id: Optional[int], guild_id: Optional[int], privileged_author_id: Optional[int]):
    bot = DiscordBot.get_instance()
    async with bot.db_pool.acquire() as conn:
        await conn.execute('''
            INSERT INTO users (discord_snowflake, coordinator_channel_ids)
            VALUES ($1, ARRAY[$2]::BIGINT[])
            ON CONFLICT (discord_snowflake) DO UPDATE
            SET coordinator_channel_ids = (
                SELECT ARRAY(SELECT DISTINCT unnest(
                    COALESCE(users.coordinator_channel_ids, ARRAY[]::BIGINT[]) || ARRAY[$2]::BIGINT[]
                ))
            ),
            updated_at = NOW()
        ''', privileged_author_id, channel_id)
    