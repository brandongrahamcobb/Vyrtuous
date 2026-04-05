from vyrtuous.bot.discord_bot import DiscordBot
import statistics

FIVE_MINUTES_SECONDS = 0.5 * 60


class SystemMonitoringService:

    __cpu_seconds = []
    __rx_bytes = []
    __tx_bytes = []

    def __init__(self):
        self.__bot = DiscordBot.get_instance()

    async def log_cpu_seconds(self):
        with open("/sys/fs/cgroup/cpu.stat", "r") as file:
            content = file.readline()
            fields = content.split()
            self.__cpu_seconds.append(int(fields[1]))

    async def log_rx_bytes(self):
        with open("/sys/class/net/eth0/statistics/rx_bytes", "r") as file:
            content = file.readline()
            self.__rx_bytes.append(int(content))

    async def log_tx_bytes(self):
        with open("/sys/class/net/eth0/statistics/tx_bytes", "r") as file:
            content = file.readline()
            self.__tx_bytes.append(int(content))

    async def calculate_cpu_usage(self):
        interval_seconds = len(self.__cpu_seconds) * FIVE_MINUTES_SECONDS
        average_difference = [
            (j - i) / 1_000_000
            for i, j in zip(self.__cpu_seconds[:-1], self.__cpu_seconds[1:])
        ]
        if len(average_difference) < 1:
            return 0.0
        usage_average = statistics.mean(average_difference)
        calculated_percentage = round((usage_average / interval_seconds) * 100, 0)
        return calculated_percentage

    async def calculate_rx_usage(self):
        interval_seconds = len(self.__rx_bytes) * FIVE_MINUTES_SECONDS
        average_difference = [
            j - i for i, j in zip(self.__rx_bytes[:-1], self.__rx_bytes[1:])
        ]
        if len(average_difference) < 1:
            return 0.0
        usage_average = statistics.mean(average_difference)
        calculated_megabytes = round(
            (((usage_average / 1024) / 1024) / interval_seconds), 0
        )
        return calculated_megabytes

    async def calculate_tx_usage(self):
        interval_seconds = len(self.__tx_bytes) * FIVE_MINUTES_SECONDS
        average_difference = [
            j - i for i, j in zip(self.__rx_bytes[:-1], self.__rx_bytes[1:])
        ]
        if len(average_difference) < 1:
            return 0.0
        usage_average = statistics.mean(average_difference)
        calculated_megabytes = round(
            (((usage_average / 1024) / 1024) / interval_seconds), 0
        )
        return calculated_megabytes
