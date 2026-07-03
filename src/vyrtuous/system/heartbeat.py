"""!/bin/python3
heartbeat.py A discord.py cog containing a heartbeat mechanism for the Vyrtuous bot.

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

from discord.ext import commands, tasks

from vyrtuous.bot.discord_bot import DiscordBot


class Heartbeat(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.heartbeat.start()

    async def cog_unload(self):
        self.heartbeat.cancel()

    @tasks.loop(seconds=5.0)
    async def heartbeat(self):
        with open("/tmp/vyrtuous_heartbeat", "w", encoding="utf-8") as f:
            if self.bot.is_ready and self.bot.latency:
                f.write("0")
            else:
                f.write("1")


async def setup(bot: DiscordBot):
    cog = Heartbeat(bot)
    await bot.add_cog(cog)
