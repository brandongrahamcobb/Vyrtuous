"""sysadmin_commands.py A discord.py cog containing sysadmin commands for the Vyrtuous bot.

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

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.bug import Bug
from vyrtuous.db.roles.developer import Developer
from vyrtuous.db.roles.sysadmin import sysadmin_predicator
from vyrtuous.fields.snowflake import (
    AppMemberSnowflake,
    MemberSnowflake,
)
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject


class SysadminCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    # DONE
    @app_commands.command(name="assign", description="Assign developer.")
    @app_commands.describe(
        reference="Include an issue reference ID",
        member="Tag a member or include their ID",
    )
    @sysadmin_predicator()
    async def assign_bug_to_developer_app_command(
        self,
        interaction: discord.Interaction,
        reference: str,
        member: AppMemberSnowflake,
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        member_dict = await do.determine_from_target(target=member)
        kwargs = member_dict.get("columns", None)

        developer = await Developer.select(**kwargs, singular=True)
        if not developer:
            return await state.end(
                warning=f"Developer not found for target ({member})."
            )

        bug = await Bug.select(id=reference, resolved=False, singular=True)
        if not bug:
            return await state.end(
                warning=f"Unresolved issue not found for reference: {reference}."
            )
        member_snowflakes = bug.member_snowflakes
        where_kwargs = {"id": bug.id}
        member_snowflakes = bug.member_snowflakes
        if developer.member_snowflake in bug.member_snowflakes:
            member_snowflakes.remove(developer.member_snowflake)
            set_kwargs = {"member_snowflakes": member_snowflakes}
            await bug.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            embed = await bug.create_embed(
                action="unassigned",
                member_snowflake=developer.member_snowflake,
                source=interaction,
            )
            return await state.end(success=embed)
        else:
            member_snowflakes.append(developer.member_snowflake)
            set_kwargs = {"member_snowflakes": member_snowflakes}
            await bug.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            embed = await bug.create_embed(
                action="assigned",
                member_snowflake=developer.member_snowflake,
                source=interaction,
            )
            await member_dict.get("object", None).send(embed=embed)
            return await state.end(success=embed)

    # DONE
    @commands.command(name="assign", help="Assign developer.")
    @sysadmin_predicator()
    async def assign_bug_to_developer_text_command(
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
        kwargs = member_dict.get("columns", None)

        developer = await Developer.select(**kwargs, singular=True)
        if not developer:
            return await state.end(
                warning=f"Developer not found for target ({member})."
            )

        bug = await Bug.select(id=reference, resolved=False, singular=True)
        if not bug:
            return await state.end(
                warning=f"Unresolved issue not found for reference: {reference}."
            )
        member_snowflakes = bug.member_snowflakes
        where_kwargs = {"id": bug.id}
        member_snowflakes = bug.member_snowflakes
        if developer.member_snowflake in bug.member_snowflakes:
            member_snowflakes.remove(developer.member_snowflake)
            set_kwargs = {"member_snowflakes": member_snowflakes}
            await bug.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            embed = await bug.create_embed(
                action="unassigned",
                member_snowflake=developer.member_snowflake,
                source=ctx,
            )
            return await state.end(success=embed)
        else:
            member_snowflakes.append(developer.member_snowflake)
            set_kwargs = {"member_snowflakes": member_snowflakes}
            await bug.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
            embed = await bug.create_embed(
                action="assigned",
                member_snowflake=developer.member_snowflake,
                source=ctx,
            )
            await member_dict.get("object", None).send(embed=embed)
            return await state.end(success=embed)

    # DONE
    @app_commands.command(name="dev", description="Grant/revoke devs.")
    @app_commands.describe(member="Tag a member or include their ID")
    @sysadmin_predicator()
    async def toggle_developer_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        action = None

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        member_dict = await do.determine_from_target(target=member)
        kwargs = member_dict.get("columns", None)

        developer = await Developer.select(**kwargs)

        if developer:
            await Developer.delete(**kwargs)
            action = "revoked"
        else:
            developer = Developer(**kwargs)
            await developer.create()
            action = "granted"

        return await state.end(
            success=f"Developer access for {member_dict.get("mention", None)} has been {action} in {interaction.guild.name}."
        )

    # DONE
    @commands.command(name="dev", help="Grant/revoke devs.")
    @sysadmin_predicator()
    async def toggle_developer_text_command(
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
        kwargs = member_dict.get("columns", None)

        developer = await Developer.select(**kwargs)

        if developer:
            await Developer.delete(**kwargs)
            action = "revoked"
        else:
            developer = Developer(**kwargs)
            await developer.create()
            action = "granted"

        return await state.end(
            success=f"Developer access for {member_dict.get("mention", None)} has been {action} in {ctx.guild.name}."
        )


async def setup(bot: DiscordBot):
    await bot.add_cog(SysadminCommands(bot))
