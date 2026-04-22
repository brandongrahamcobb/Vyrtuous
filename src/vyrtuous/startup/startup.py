# startup.py Upon loading the bot, in memory data structures need to be populated for quicker memory access instead of running a DB query each time.
from discord.ext import commands
from vyrtuous.active_members.active_member_service import ActiveMemberService
from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bug.bug_service import BugService
from vyrtuous.duration.duration_builder import DurationBuilder
from vyrtuous.ban.ban_service import BanService
from vyrtuous.flag.flag_service import FlagService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.data_service import DataService
from vyrtuous.utils.message_service import PaginatorService
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService
from vyrtuous.vegan.vegan_service import VeganService


class Startup(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__author_service = AuthorService()
        self.__bot = bot
        self.__emoji = Emojis()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__active_member_service = ActiveMemberService(
            bot=self.__bot, database_factory=self.__database_factory
        )
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__duration_builder = DurationBuilder()
        self.__paginator_service = PaginatorService(bot=self.__bot)
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__administrator_service = AdministratorService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__coordinator_service = CoordinatorService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__developer_service = DeveloperService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            database_factory=self.__database_factory,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
        )
        self.__guild_owner_service = GuildOwnerService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__sysadmin_service = SysadminService(
            active_member_service=self.__active_member_service,
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__moderator_service = ModeratorService(
            active_member_service=self.__active_member_service,
            administrator_service=self.__administrator_service,
            author_service=self.__author_service,
            bot=self.__bot,
            coordinator_service=self.__coordinator_service,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            guild_owner_service=self.__guild_owner_service,
            sysadmin_service=self.__sysadmin_service,
        )
        self.__data_service = DataService(
            database_factory=self.__database_factory,
            duration_builder=self.__duration_builder,
            moderator_service=self.__moderator_service,
        )
        self.__stream_service = StreamService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            paginator_service=self.__paginator_service,
        )
        self.__voice_mute_service = VoiceMuteService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            stream_service=self.__stream_service,
        )
        self.__ban_service = BanService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
        )
        self.__flag_service = FlagService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__vegan_service = VeganService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__text_mute_service = TextMuteService(
            active_member_service=self.__active_member_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            data_service=self.__data_service,
            dictionary_service=self.__dictionary_service,
            duration_builder=self.__duration_builder,
            emoji=self.__emoji,
        )

    async def cog_load(self):
        await self.__moderator_service.populate()
        await self.__coordinator_service.populate()
        await self.__administrator_service.populate()
        await self.__guild_owner_service.populate()
        await self.__developer_service.populate()
        await self.__sysadmin_service.populate()
        await self.__active_member_service.populate()
        await self.__ban_service.populate()
        await self.__flag_service.populate()
        await self.__text_mute_service.populate()
        await self.__vegan_service.populate()
        await self.__voice_mute_service.populate()


async def setup(bot: DiscordBot):
    await bot.add_cog(Startup(bot=bot))
