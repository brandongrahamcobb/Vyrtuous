"""system_owner_commands.py A discord.py cog containing system owner commands for the Vyrtuous bot.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.logs.developer_log import DeveloperLog
from vyrtuous.database.roles.developer import Developer
from vyrtuous.properties.snowflake import (
    AppMemberSnowflake,
    MemberSnowflake,
)
from vyrtuous.service.check_service import (
    sys_owner_predicator,
)
from vyrtuous.service.logging_service import logger
from vyrtuous.service.messaging.message_service import MessageService
from vyrtuous.service.messaging.state_service import StateService
from vyrtuous.service.resolution.discord_object_service import DiscordObject
from vyrtuous.utils.emojis import get_random_emoji


class SystemOwnerCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    # DONE
    @app_commands.command(name="assign", description="Assign developer.")
    @app_commands.describe(
        reference="Include an issue reference ID",
        member="Tag a member or include their ID",
    )
    @sys_owner_predicator()
    async def toggle_issue_to_developer_app_command(
        self,
        interaction: discord.Interaction,
        reference: str,
        member: AppMemberSnowflake,
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        member_dict = await do.determine_from_target(target=member)
        kwargs = member_dict["columns"]

        developer = await Developer.select(**kwargs)

        developer_log = await DeveloperLog.select(id=reference, resolved=False)
        if developer_log:
            channel = self.bot.get_channel(developer_log.channel_snowflake)
            try:
                msg = await channel.fetch_message(developer_log.message_snowflake)
                link = msg.jump_url
            except discord.NotFound:
                try:
                    return await state.end(
                        warning=f"Message reference not found: {reference}."
                    )
                except Exception as e:
                    link = "Unknown message"
                    logger.warning(str(e).capitalize())
            if developer.member_snowflake in developer_log.developer_snowflakes:
                await developer_log.unassign(member_snowflake=member_dict["id"])
                return await state.end(
                    success=f"Developer {member_dict['mention']} unassigned for issue by {interaction.user.mention}: {link}\n**Notes:** {developer_log.notes}."
                )
            else:
                await developer_log.assign(member_snowflake=member_dict["id"])
                await state.end(
                    success=f"Developer {member_dict['mention']} assigned to issue by {interaction.user.mention}: {link}\n**Notes:** {developer_log.notes}."
                )
                return await member_dict["object"].send(
                    f"{get_random_emoji()} Developer {member_dict['mention']} assigned to issue by {interaction.user.mention}: {link}\n**Notes:** {developer_log.notes}"
                )
        else:
            return await state.end(
                warning=f"Unresolved issue not found for reference: {reference}."
            )

    # DONE
    @commands.command(name="assign", help="Assign developer.")
    @sys_owner_predicator()
    async def toggle_issue_to_developer_text_command(
        self,
        ctx: commands.Context,
        reference: str = commands.parameter(
            description="Include an issue reference ID"
        ),
        member: MemberSnowflake = commands.parameter(
            description="Tag a member or include their ID"
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        member_dict = await do.determine_from_target(target=member)
        kwargs = member_dict["columns"]

        developer = await Developer.select(**kwargs)

        developer_log = await DeveloperLog.select(id=reference, resolved=False)
        if developer_log:
            channel = self.bot.get_channel(developer_log.channel_snowflake)
            try:
                msg = await channel.fetch_message(developer_log.message_snowflake)
                link = msg.jump_url
            except discord.NotFound:
                try:
                    return await state.end(
                        warning=f"Message reference not found: {reference}."
                    )
                except Exception as e:
                    link = "Unknown message"
                    logger.warning(str(e).capitalize())
            if developer.member_snowflake in developer_log.developer_snowflakes:
                await developer_log.unassign(member_snowflake=member_dict["id"])
                return await state.end(
                    success=f"Developer {member_dict['mention']} unassigned for issue by {ctx.author.mention}: {link}\n**Notes:** {developer_log.notes}."
                )
            else:
                await developer_log.assign(member_snowflake=member_dict["id"])
                await state.end(
                    success=f"Developer {member_dict['mention']} assigned for issue by {ctx.author.mention}: {link}\n**Notes:** {developer_log.notes}."
                )
                return await member_dict["object"].send(
                    f"{get_random_emoji()} Developer {member_dict['mention']} assigned for issue by {ctx.author.mention}: {link}\n**Notes:** {developer_log.notes}"
                )
        else:
            return await state.end(
                warning=f"Unresolved issue not found for reference: {reference}."
            )

    # DONE
    @app_commands.command(name="dev", description="Grant/revoke devs.")
    @app_commands.describe(member="Tag a member or include their ID")
    @sys_owner_predicator()
    async def create_developer_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        action = None

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        member_dict = await do.determine_from_target(target=member)
        kwargs = member_dict["columns"]

        developer = await Developer.select(**kwargs)

        if developer:
            await Developer.delete(**kwargs)
            action = "revoked"
        else:
            developer = Developer(**kwargs)
            await developer.create()
            action = "granted"

        return await state.end(
            success=f"Developer access for {member_dict['mention']} has been {action} in {interaction.guild.name}."
        )

    # DONE
    @commands.command(name="dev", help="Grant/revoke devs.")
    @sys_owner_predicator()
    async def create_developer_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Tag a member or include their ID"
        ),
    ):
        action = None

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)

        member_dict = await do.determine_from_target(target=member)
        kwargs = member_dict["columns"]

        developer = await Developer.select(**kwargs)

        if developer:
            await Developer.delete(**kwargs)
            action = "revoked"
        else:
            developer = Developer(**kwargs)
            await developer.create()
            action = "granted"

        return await state.end(
            success=f"Developer access for {member_dict['mention']} has been {action} in {ctx.guild.name}."
        )


async def setup(bot: DiscordBot):
    await bot.add_cog(SystemOwnerCommands(bot))
