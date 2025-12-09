from discord.ext import tasks, commands

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
            if (self.bot.is_ready and self.bot.latency):
                f.write("0")
            else:
                f.write("1")

async def setup(bot: DiscordBot):
    cog = Heartbeat(bot)
    await bot.add_cog(cog)
