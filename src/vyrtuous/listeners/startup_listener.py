"""!/bin/python3
startup_listener.py The purpose of this program is to provide a cog when loading the bot, in memory data structures need to be populated for quicker memory access instead of running a DB query each time.

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

from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.utils.channels import video_channel_service
from vyrtuous.utils.moderation import flag_service, voice_mute_service
from vyrtuous.utils.permissions import permission_service
from vyrtuous.utils.users import active_member_service


class Startup(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def cog_load(self) -> None:
        await active_member_service.populate()
        await flag_service.populate()
        permission_service.populate(
            permission_state=self.__bot.registry.get(PermissionState)
        )
        await video_channel_service.populate()
        await voice_mute_service.populate(target="auto")


async def setup(bot: DiscordBot):
    await bot.add_cog(Startup(bot))
