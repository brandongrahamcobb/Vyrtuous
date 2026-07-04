"""!/bin/python3
dev_text_commands.py A discord.py cog containing developer commands for the Vyrtuous bot.

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

from datetime import timedelta

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.developer import NotDeveloper
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.statistics import system_monitoring_service
from vyrtuous.utils.users import developer_service, sysadmin_service


class DeveloperTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Developer"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx: commands.Context) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        context = SnowflakeContext(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        for verify in (
            sysadmin_service.is_sysadmin_wrapper,
            developer_service.is_developer_wrapper,
        ):
            try:
                if await verify(context=context):
                    return True
            except commands.CheckFailure:
                continue
        raise NotDeveloper

    @commands.command(name="ping", help="Ping me!")
    async def ping_text_command(self, ctx: commands.Context) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        return await tick.end(success="Pong!")

    @commands.command(name="stats", help="Lists stats.")
    async def list_stats_text_command(self, ctx: commands.Context) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        embed = discord.Embed(title="Statistics")
        cpu_usage = await system_monitoring_service.calculate_cpu_usage()
        embed.add_field(name="CPU %", value=cpu_usage, inline=True)
        with open("/sys/fs/cgroup/memory.current", "r") as file:
            bits = int(file.read())
            memory_usage = round((bits / 1024) / 1024, 0)
            embed.add_field(name="RAM", value=f"{memory_usage} MB", inline=True)
        with open("/proc/uptime", "r") as file:
            content = file.readline()
            fields = content.split()
            uptime_seconds = float(fields[0])
            time = timedelta(seconds=uptime_seconds)
        embed.add_field(name="Uptime", value=f"{str(time)}", inline=True)
        rx_usage = await system_monitoring_service.calculate_rx_usage()
        embed.add_field(name="RX MB", value=f"{rx_usage} MB", inline=True)
        tx_usage = await system_monitoring_service.calculate_tx_usage()
        embed.add_field(name="TX MB", value=f"{tx_usage} MB", inline=True)
        number_of_servers = len(self.__bot.guilds)
        embed.add_field(name="Servers", value=number_of_servers, inline=True)
        return await tick.end(success=embed)


async def setup(bot: DiscordBot):
    await bot.add_cog(DeveloperTextCommands(bot=bot))
