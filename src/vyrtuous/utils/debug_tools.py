# from vyrtuous.bot.discord_bot import DiscordBot
#
#
# class DebugTools(commands.Cog)
#
#     def __init__(self):
#         self.__bot = DiscordBot()
#
#     @commands.command(name="")
#
#
# async def setup(bot):
#     await bot.add_cog(DebugTools(bot))
import discord
from vyrtuous.active_members.active_member_service import ActiveMemberService
from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_builder import DurationBuilder
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.data_service import DataService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.emojis import Emojis


class DebugTools:

    def __init__(self, bot: DiscordBot | None):
        self.__bot = bot
        self.__author_service = AuthorService()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__active_member_service = ActiveMemberService(
            bot=self.__bot, database_factory=self.__database_factory
        )
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__duration_builder = DurationBuilder()
        # self.__ban_service = BanService(
        #     active_member_service=self.__active_member_service,
        #     bot=self.__bot,
        #     database_factory=self.__database_factory,
        #     data_service=self.__data_service,
        #     dictionary_service=self.__dictionary_service,
        #     duration_builder=self.__duration_builder,
        #     emoji=self.__emoji,
        # )
        # self.__text_mute_service = TextMuteService(
        #     active_member_service=self.__active_member_service,
        #     bot=self.__bot,
        #     database_factory=self.__database_factory,
        #     data_service=self.__data_service,
        #     dictionary_service=self.__dictionary_service,
        #     duration_builder=self.__duration_builder,
        #     emoji=self.__emoji,
        # )
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

    async def cleanup_stale_overwrites(self, guild: discord.Guild):
        for channel in guild.channels:
            for obj, overwrite in channel.overwrites.items():
                if isinstance(obj, discord.Role):
                    continue
                elif isinstance(obj, discord.Member):
                    continue
