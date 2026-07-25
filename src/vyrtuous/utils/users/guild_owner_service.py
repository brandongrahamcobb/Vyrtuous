"""!/bin/python3
guild_owner_service.py The purpose of this program is to service the guild owners.

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
from vyrtuous.cache.guild_owner import GuildOwner, NotGuildOwner

MODEL = GuildOwner


async def is_guild_owner(
    member_snowflake: int, guild_snowflake: int | None = None
) -> bool:
    bot: DiscordBot = DiscordBot.get_instance()
    if guild_snowflake is None:
        for guild in bot.guilds:
            if guild and guild.owner_id == member_snowflake:
                return True
    else:
        guild = bot.get_guild(guild_snowflake)
        if guild and guild.owner_id == member_snowflake:
            return True
    raise NotGuildOwner(guild_snowflake=guild_snowflake)
