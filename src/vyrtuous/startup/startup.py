# startup.py Upon loading the bot, in memory data structures need to be populated for quicker memory access instead of running a DB query each time.
from discord.ext import commands
from vyrtuous.active_members.active_member_service import ActiveMemberService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService


class Startup(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__active_member_service = ActiveMemberService(
            bot=self.__bot, database_factory=self.__database_factory
        )

    async def cog_load(self):
        await self.__active_member_service.populate()


async def setup(bot: DiscordBot):
    await bot.add_cog(Startup(bot=bot))
