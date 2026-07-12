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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import (
    list_aliases,
    list_bans,
    list_coordinators,
    list_flags,
    list_moderators,
    list_text_mutes,
    list_voice_mutes,
)
from vyrtuous.models import multi_converter
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import moderator_service
from vyrtuous.view.flag_view import FlagView
from vyrtuous.view.infraction_view import InfractionView
from vyrtuous.view.modify_infraction_view import ModifyInfractionView
from vyrtuous.view.view_context import ViewContext


class ModeratorAppCommands(commands.Cog):

    PERMISSION_LEVEL = "Moderator"

    def __init__(
        self,
        bot: DiscordBot,
    ):
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

    @app_commands.command(name="duration", description="Modify a duration.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def change_moderation_duration_app_command(
        self, interaction: discord.Interaction, member: str
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        ctx = ViewContext(interaction=interaction, member_snowflake=member_snowflake)
        view = ModifyInfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            modal="duration",
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a category", view=view, ephemeral=True
        )

    @app_commands.command(name="reason", description="Modify a reason.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def change_moderation_reason_app_command(
        self, interaction: discord.Interaction, member: str
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        ctx = ViewContext(interaction=interaction, member_snowflake=member_snowflake)
        view = ModifyInfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            modal="reason",
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a category", view=view, ephemeral=True
        )

    @app_commands.command(name="ban", description="Create a ban.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_ban_app_command(
        self, interaction: discord.Interaction, member: str
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        ctx = ViewContext(interaction=interaction, member_snowflake=member_snowflake)
        ctx.category = "ban"
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )

    @app_commands.command(name="flag", description="Create a flag.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_flag_app_command(
        self, interaction: discord.Interaction, member: str
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        ctx = ViewContext(interaction=interaction, member_snowflake=member_snowflake)
        ctx.category = "flag"
        view = FlagView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel", view=view, ephemeral=True
        )

    @app_commands.command(name="mute", description="Create a voice mute.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_voice_mute_app_command(
        self, interaction: discord.Interaction, member: str
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        ctx = ViewContext(interaction=interaction, member_snowflake=member_snowflake)
        ctx.category = "vmute"
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )

    @app_commands.command(name="tmute", description="Create a text-mute.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def create_app_mute_app_command(
        self, interaction: discord.Interaction, member: str
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        ctx = ViewContext(interaction=interaction, member_snowflake=member_snowflake)
        ctx.category = "tmute"
        view = InfractionView(
            author_snowflake=interaction.user.id,
            ctx=ctx,
            tick=tick,
        )
        await view.setup()
        await interaction.response.send_message(
            content="Select a channel and a duration", view=view, ephemeral=True
        )

    @app_commands.command(name="bans", description="List bans.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_bans_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_bans.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="cmds", description="List aliases.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, or server ID."
    )
    async def list_commands_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_aliases.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="coords", description="Lists coords.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_coordinators_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_coordinators.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="flags", description="List flags.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_flags_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_flags.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="mods", description="Lists mods.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_moderators_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_moderators.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="mutes", description="List mutes.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_voice_mutes.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="summary", description="List user moderation.")
    @app_commands.describe(
        member="Specify a member ID/mention.",
    )
    async def list_moderation_summary_app_command(
        self, interaction: discord.Interaction, member: str
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        obj = multi_converter.transform(interaction=interaction, argument=member)
        pages: list[discord.Embed] = []
        is_at_home = at_home(source=interaction)
        services = []
        services.append(ban_service)
        services.append(flag_service)
        services.append(text_mute_service)
        services.append(voice_mute_service)
        for service in services:
            summary_pages = await service.build_pages(obj=obj, is_at_home=is_at_home)
            if isinstance(summary_pages, list):
                for page in summary_pages:
                    if isinstance(page, discord.Embed):
                        pages.append(page)
        if not pages:
            return await tick.end(success="No infractions found")
        return await tick.end(success=pages)

    @app_commands.command(name="tmutes", description="List text-mutes.")
    @app_commands.describe(
        target="Specify one of: `all`, a channel ID/mention, member ID/mention or server ID.",
    )
    async def list_text_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_text_mutes.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(ModeratorAppCommands(bot=bot))
