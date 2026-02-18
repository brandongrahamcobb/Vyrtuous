"""!/bin/python3
moderator_app_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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

from typing import Any, Coroutine, Union

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.administrator.administrator_service import AdministratorService
from vyrtuous.ban.ban import Ban
from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.cap.cap_service import CapService
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.developer.developer_service import DeveloperService
from vyrtuous.duration.duration_service import DurationService
from vyrtuous.flag.flag_service import FlagService
from vyrtuous.modal.duration_modal import DurationModal
from vyrtuous.modal.reason_modal import ReasonModal
from vyrtuous.moderator.moderator import Moderator
from vyrtuous.moderator.moderator_service import ModeratorService, NotModerator
from vyrtuous.owner.guild_owner_service import GuildOwnerService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.text_mute.text_mute import TextMute
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import DiscordObjectService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.state_service import StateService
from vyrtuous.utils.view_context import ViewContext
from vyrtuous.view.data_view import DataView
from vyrtuous.view.infraction_view import InfractionView
from vyrtuous.view.modify_infraction_view import ModifyInfractionView
from vyrtuous.voice_mute.voice_mute import VoiceMute
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService


class ModeratorAppCommands(commands.Cog):
    ROLE = Moderator

    def __init__(
        self,
        *,
        bot: DiscordBot | None = None,
    ):
        self.__author_service = AuthorService()
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__emoji = Emojis()
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__duration_service = DurationService()
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__developer_service = DeveloperService(
            author_service=self.__author_service,
            bug_service=self.__bug_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            emoji=self.__emoji,
        )
        self.__moderator_service = ModeratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__coordinator_service = CoordinatorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__administrator_service = AdministratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__guild_owner_service = GuildOwnerService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__stream_service = StreamService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
        )
        self.__voice_mute_service = VoiceMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            stream_service=self.__stream_service,
        )
        self.__ban_service = BanService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
        )
        self.__flag_service = FlagService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__text_mute_service = TextMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__cap_service = CapService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
        )
        self.__discord_object_service = DiscordObjectService()

    async def cog_check(self, interaction) -> Coroutine[Any, Any, bool]:
        async def predicate(
            source: Union[commands.Context, discord.Interaction, discord.Message],
        ):
            for verify in (
                self.__sysadmin_service.is_sysadmin_wrapper,
                self.__developer_service.is_developer_wrapper,
                self.__guild_owner_service.is_guild_owner_wrapper,
                self.__administrator_service.is_administrator_wrapper,
                self.__coordinator_service.is_coordinator_at_all_wrapper,
                self.__moderator_service.is_moderator_at_all_wrapper,
            ):
                try:
                    if await verify(source):
                        return True
                except app_commands.CheckFailure:
                    continue
            raise NotModerator

        predicate._permission_level = "Moderator"
        return await predicate(interaction)

    @app_commands.command(name="data", description="Create a chart.")
    async def create_data_app_command(self, interaction: discord.Interaction):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            interaction=interaction,
        )
        default_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        view = DataView(
            ban_service=self.__ban_service,
            duration_service=self.__duration_service,
            flag_service=self.__flag_service,
            interaction=interaction,
            moderator_service=self.__moderator_service,
            text_mute_service=self.__text_mute_service,
            voice_mute_service=self.__voice_mute_service,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel, duration and infraction",
            view=view,
            ephemeral=True,
        )

    @app_commands.command(name="duration", description="Modify a duration.")
    @app_commands.describe(member="The ID or mention of the member.")
    async def change_moderation_duration_app_command(
        self,
        interaction: discord.Interaction,
        member: discord.Member = app_commands.parameter(
            converter=app_commands.MemberConverter,
            default=None,
            description="Tag a member or include their ID",
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            interaction=interaction,
        )
        member_dict = self.__discord_object_service.to_dict(obj=member)
        ctx = ViewContext(
            ban_service=self.__ban_service,
            flag_service=self.__flag_service,
            interaction=interaction,
            moderator_service=self.__moderator_service,
            text_mute_service=self.__text_mute_service,
            voice_mute_service=self.__voice_mute_service,
        )
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = ModifyInfractionView(
            ctx=ctx,
            modal=DurationModal,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and an infraction", view=view, ephemeral=True
        )

    @app_commands.command(name="reason", description="Modify a reason.")
    @app_commands.describe(member="The ID or mention of the member.")
    async def change_moderation_reason_app_command(
        self,
        interaction: discord.Interaction,
        member: discord.Member = app_commands.parameter(
            converter=app_commands.MemberConverter,
            default=None,
            description="Tag a member or include their ID",
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            interaction=interaction,
        )
        default_kwargs = {
            "channel_snowflake": int(interaction.channel.id),
            "guild_snowflake": int(interaction.guild.id),
            "member_snowflake": int(interaction.user.id),
        }
        member_dict = self.__discord_object_service.to_dict(obj=member)
        ctx = ViewContext(
            ban_service=self.__ban_service,
            flag_service=self.__flag_service,
            interaction=interaction,
            moderator_service=self.__moderator_service,
            text_mute_service=self.__text_mute_service,
            voice_mute_service=self.__voice_mute_service,
        )
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = ModifyInfractionView(
            ctx=ctx,
            modal=ReasonModal,
            state=state,
            ban_service=self.__ban_service,
            flag_service=self.__flag_service,
            text_mute_service=self.__text_mute_service,
            voice_mute_service=self.__voice_mute_service,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and an infraction", view=view, ephemeral=True
        )

    @app_commands.command(name="vban", description="Create a ban.")
    @app_commands.describe(member="The ID or mention of the member.")
    async def create_ban_app_command(
        self,
        interaction: discord.Interaction,
        member: discord.Member = app_commands.parameter(
            converter=app_commands.MemberConverter,
            default=None,
            description="Tag a member or include their ID",
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            interaction=interaction,
        )
        member_dict = self.__discord_object_service.to_dict(obj=member)
        ctx = ViewContext(
            ban_service=self.__ban_service,
            flag_service=self.__flag_service,
            interaction=interaction,
            moderator_service=self.__moderator_service,
            text_mute_service=self.__text_mute_service,
            voice_mute_service=self.__voice_mute_service,
        )
        ctx.infraction = Ban
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = InfractionView(
            cap_service=self.__cap_service,
            ctx=ctx,
            modal=ReasonModal,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )

    @app_commands.command(name="vmute", description="Create a mute.")
    @app_commands.describe(member="The ID or mention of the member.")
    async def create_voice_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: discord.Member = app_commands.parameter(
            converter=app_commands.MemberConverter,
            default=None,
            description="Tag a member or include their ID",
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            interaction=interaction,
        )
        member_dict = self.__discord_object_service.to_dict(obj=member)
        ctx = ViewContext(
            ban_service=self.__ban_service,
            flag_service=self.__flag_service,
            interaction=interaction,
            moderator_service=self.__moderator_service,
            text_mute_service=self.__text_mute_service,
            voice_mute_service=self.__voice_mute_service,
        )
        ctx.infraction = VoiceMute
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = InfractionView(
            cap_service=self.__cap_service,
            ctx=ctx,
            modal=ReasonModal,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )

    @app_commands.command(name="vtmute", description="Create a text-mute.")
    @app_commands.describe(member="The ID or mention of the member.")
    async def create_text_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: discord.Member = app_commands.parameter(
            converter=app_commands.MemberConverter,
            default=None,
            description="Tag a member or include their ID",
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
            interaction=interaction,
        )
        member_dict = self.__discord_object_service.to_dict(obj=member)
        ctx = ViewContext(
            ban_service=self.__ban_service,
            flag_service=self.__flag_service,
            interaction=interaction,
            moderator_service=self.__moderator_service,
            text_mute_service=self.__text_mute_service,
            voice_mute_service=self.__voice_mute_service,
        )
        ctx.infraction = TextMute
        await ctx.setup(target_member_snowflake=member_dict.get("id", None))
        view = InfractionView(
            cap_service=self.__cap_service,
            ctx=ctx,
            modal=ReasonModal,
            state=state,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )


async def setup(bot: DiscordBot):
    cog = ModeratorAppCommands(bot)
    await bot.add_cog(cog)
