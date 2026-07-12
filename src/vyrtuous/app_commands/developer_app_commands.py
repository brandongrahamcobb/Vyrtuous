"""!/bin/python3
dev_app_commands.py A discord.py cog containing developer commands for the Vyrtuous bot.

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
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.statistics import system_monitoring_service
from vyrtuous.utils.users import moderator_service


class DeveloperAppCommands(commands.Cog):

    PERMISSION_LEVEL = "Developer"

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def interaction_check(self, interaction: discord.Interaction):
        if interaction.guild is None:
            raise commands.CheckFailure("This command must be used inside a server.")
        if interaction.channel is None:
            raise commands.CheckFailure("This command must be used in a valid channel.")
        await moderator_service.check_minimum_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @app_commands.command(name="ping", description="Ping me!")
    async def ping_app_command(
        self, interaction: discord.Interaction
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        return await tick.end(success="Pong!")

    @app_commands.command(name="stats", description="Lists stats.")
    async def list_stats_app_command(
        self, interaction: discord.Interaction
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
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
    await bot.add_cog(DeveloperAppCommands(bot=bot))
