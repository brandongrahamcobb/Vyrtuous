''' cap.py The purpose of this program is to provide the Cap utility class.
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

class Cap:
 
    @classmethod
    async def fetch_caps_for_channel(self, guild_id: int, channel_id: int) -> list[tuple[str, str]]:
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                'SELECT duration_seconds, moderation_type FROM active_caps WHERE guild_id=$1 AND channel_id=$2',
                guild_id, channel_id
            )
            return [(r['duration_seconds'], r['moderation_type']) for r in rows]

    @classmethod
    async def delete_cap_by_channel_moderation_type_and_room_name(self, channel_id: Optional[int], guild_id: Optional[int], moderation_type: Optional[str], room_name: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
              await conn.execute('''
                DELETE FROM active_caps
                WHERE channel_id=$1 AND guild_id=$2 AND moderation_type=$3 AND room_name=$4
            ''', channel_id, guild_id, moderation_type, room_name)