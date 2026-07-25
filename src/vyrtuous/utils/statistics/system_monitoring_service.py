"""!/bin/python3
system_monitoring_service.py The purpose of this program is to service system statistics.

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

import statistics

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import SystemResourcesState

FIVE_MINUTES_SECONDS = 5 * 60


async def log_cpu_seconds() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    with open("/sys/fs/cgroup/cpu.stat", "r") as file:
        content = file.readline()
        fields = content.split()
        bot.registry.get(SystemResourcesState).cpu_seconds.append(int(fields[1]))


async def log_rx_bytes() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    with open("/sys/class/net/eth0/statistics/rx_bytes", "r") as file:
        content = file.readline()
        bot.registry.get(SystemResourcesState).rx_bytes.append(int(content))


async def log_tx_bytes() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    with open("/sys/class/net/eth0/statistics/tx_bytes", "r") as file:
        content = file.readline()
        bot.registry.get(SystemResourcesState).tx_bytes.append(int(content))


async def calculate_cpu_usage() -> float:
    bot: DiscordBot = DiscordBot.get_instance()
    state = bot.registry.get(SystemResourcesState)
    interval_seconds = len(state.cpu_seconds) * FIVE_MINUTES_SECONDS
    average_difference = [
        (j - i) / 1_000_000
        for i, j in zip(state.cpu_seconds[:-1], state.cpu_seconds[1:])
    ]
    if len(average_difference) < 1:
        return 0.0
    usage_average = statistics.mean(average_difference)
    calculated_percentage = round((usage_average / interval_seconds) * 100, 0)
    return calculated_percentage


async def calculate_rx_usage() -> float:
    bot: DiscordBot = DiscordBot.get_instance()
    state = bot.registry.get(SystemResourcesState)
    interval_seconds = len(state.rx_bytes) * FIVE_MINUTES_SECONDS
    average_difference = [
        j - i for i, j in zip(state.rx_bytes[:-1], state.rx_bytes[1:])
    ]
    if len(average_difference) < 1:
        return 0.0
    usage_average = statistics.mean(average_difference)
    calculated_megabytes = round(
        (((usage_average / 1024) / 1024) / interval_seconds), 0
    )
    return calculated_megabytes


async def calculate_tx_usage() -> float:
    bot: DiscordBot = DiscordBot.get_instance()
    state = bot.registry.get(SystemResourcesState)
    interval_seconds = len(state.tx_bytes) * FIVE_MINUTES_SECONDS
    average_difference = [
        j - i for i, j in zip(state.rx_bytes[:-1], state.rx_bytes[1:])
    ]
    if len(average_difference) < 1:
        return 0.0
    usage_average = statistics.mean(average_difference)
    calculated_megabytes = round(
        (((usage_average / 1024) / 1024) / interval_seconds), 0
    )
    return calculated_megabytes
