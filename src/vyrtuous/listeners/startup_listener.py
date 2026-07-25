# startup.py Upon loading the bot, in memory data structures need to be populated for quicker memory access instead of running a DB query each time.
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
