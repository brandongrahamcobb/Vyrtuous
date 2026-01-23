"""guild_owner_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

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
from vyrtuous.db.roles.administrator import Administrator, AdministratorRole
from vyrtuous.db.roles.guild_owner import guild_owner_predicator
from vyrtuous.fields.snowflake import (
    AppMemberSnowflake,
    AppRoleSnowflake,
    MemberSnowflake,
    RoleSnowflake,
)
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import DiscordObject

from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.invincibility import Invincibility


class GuildOwnerCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    # DONE
    @app_commands.command(name="arole", description="Role -> Administrator.")
    @guild_owner_predicator()
    async def grant_administrator_by_role_app_command(
        self, interaction: discord.Interaction, role: AppRoleSnowflake
    ):
        action = None
        chunk_size, field_count, pages = 7, 0, []
        title = f"{get_random_emoji()} Administrators and Roles"

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)

        role_dict = await do.determine_from_target(target=role)
        kwargs = role_dict.get("columns", None)

        administrators = await Administrator.select(
            role_snowflakes=role_dict.get("id", None), inside=True
        )
        administrator_roles = await AdministratorRole.select(
            role_snowflake=role_dict.get("id", None)
        )

        if administrator_roles:
            for administrator_role in administrator_roles:
                await AdministratorRole.delete(**kwargs)
            action = "revoked"
        else:
            action = "granted"
        if administrators:
            for administrator in administrators:
                revoked_members = {}
                member = interaction.guild.get_member(administrator.member_snowflake)
                await Administrator.delete(
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=administrator.member_snowflake,
                )
                revoked_members.setdefault(administrator.guild_snowflake, {}).setdefault(
                    role_dict.get("id", None), []
                ).append(member)
        else:
            revoked_members = {}

        if action != "revoked":
            granted_members = {}
            granted_members.setdefault(interaction.guild.id, {})[
                role_dict.get("id", None)
            ] = []
            administrator_role = AdministratorRole(**kwargs)
            await administrator_role.create()
            role_snowflakes = [role_dict.get("id", None)]
            for member in role_dict.get("object", None).members:
                administrator = Administrator(
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=member.id,
                    role_snowflakes=role_snowflakes,
                )
                await administrator.create()
                granted_members[interaction.guild.id][role_dict.get("id", None)].append(
                    member
                )

        embed = discord.Embed(
            title=title,
            description=f"`{role_dict.get('name', None)}` was `{action}`.",
            color=discord.Color.red() if action == "revoked" else discord.Color.green(),
        )
        embed.add_field(
            name="Role ID", value=str(role_dict.get("id", None)), inline=False
        )
        embed.add_field(name="Guild", value=interaction.guild.name, inline=False)
        pages.append(embed)
        if action == "revoked":
            members = revoked_members.get(interaction.guild.id, {}).get(
                role_dict.get("id", None), []
            )
        else:
            members = granted_members.get(interaction.guild.id, {}).get(
            role_dict.get("id", None), []
        )
    
        chunks = []
        chunk = []
        for member in members:
            chunk.append(member)
            if len(chunk) == chunk_size:
                chunks.append(chunk)
                chunk = []
        if chunk:
            chunks.append(chunk)
        field_count = 1
        page_number = 1
        for chunk in chunks:
            embed = discord.Embed(
                title=f"Members {action.capitalize()}",
                color=(
                    discord.Color.red()
                    if action == "revoked"
                    else discord.Color.green()
                ),
            )
            for member in chunk:
                embed.add_field(
                    name=f"{field_count}. {member}",
                    value=f"{member.mention} ({member.id})",
                    inline=False,
                )
                field_count += 1
            embed.set_footer(text=f"Page {page_number}")
            pages.append(embed)
            page_number += 1

        await StateService.send_pages(obj=AdministratorRole, pages=pages, state=state)

    # DONE
    @commands.command(name="arole", help="Role -> Administrator.")
    @guild_owner_predicator()
    async def grant_administrator_by_role_text_command(
        self, ctx: commands.Context, role: RoleSnowflake
    ):
        action = None
        chunk_size, field_count, pages = 7, 0, []
        title = f"{get_random_emoji()} Administrators and Roles"
    
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
    
        role_dict = await do.determine_from_target(target=role)
        kwargs = role_dict.get("columns", None)
    
        administrators = await Administrator.select(
            role_snowflakes=role_dict.get("id", None), inside=True
        )
        administrator_roles = await AdministratorRole.select(
            role_snowflake=role_dict.get("id", None)
        )
        if administrator_roles:
            for administrator_role in administrator_roles:
                await AdministratorRole.delete(**kwargs)
            action = "revoked"
        else:
            action = "granted"
        if administrators:
            for administrator in administrators:
                revoked_members = {}
                member = ctx.guild.get_member(administrator.member_snowflake)
                await Administrator.delete(
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=administrator.member_snowflake,
                )
                revoked_members.setdefault(administrator.guild_snowflake, {}).setdefault(
                    role_dict.get("id", None), []
                ).append(member)
        else:
            revoked_members = {}
        if action != "revoked":
            granted_members = {}
            granted_members.setdefault(ctx.guild.id, {})[role_dict.get("id", None)] = []
            administrator_role = AdministratorRole(**kwargs)
            await administrator_role.create()
            role_snowflakes = [role_dict.get("id", None)]
            for member in role_dict.get("object", None).members:
                administrator = Administrator(
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=member.id,
                    role_snowflakes=role_snowflakes,
                )
                await administrator.create()
                granted_members[ctx.guild.id][role_dict.get("id", None)].append(member)
    
        embed = discord.Embed(
            title=title,
            description=f"`{role_dict.get('name', None)}` was `{action}`.",
            color=discord.Color.red() if action == "revoked" else discord.Color.green(),
        )
        embed.add_field(
            name="Role ID", value=str(role_dict.get("id", None)), inline=False
        )
        embed.add_field(name="Guild", value=ctx.guild.name, inline=False)
        pages.append(embed)
        if action == "revoked":
            members = revoked_members.get(ctx.guild.id, {}).get(
                role_dict.get("id", None), []
            )
        else:
            members = granted_members.get(ctx.guild.id, {}).get(
            role_dict.get("id", None), []
        )
        chunks = []
        chunk = []
        for member in members:
            chunk.append(member)
            if len(chunk) == chunk_size:
                chunks.append(chunk)
                chunk = []
        if chunk:
            chunks.append(chunk)
        field_count = 1
        page_number = 1
        for chunk in chunks:
            embed = discord.Embed(
                title=f"Members {action.capitalize()}",
                color=(
                    discord.Color.red()
                    if action == "revoked"
                    else discord.Color.green()
                ),
            )
            for member in chunk:
                embed.add_field(
                    name=f"{field_count}. {member}",
                    value=f"{member.mention} ({member.id})",
                    inline=False,
                )
                field_count += 1
            embed.set_footer(text=f"Page {page_number}")
            pages.append(embed)
            page_number += 1
    
        await StateService.send_pages(obj=AdministratorRole, pages=pages, state=state)

    # DONE
    @app_commands.command(name="hero", description="Grant/revoke invincibility.")
    @app_commands.describe(member="Tag a member or include their ID")
    @guild_owner_predicator()
    async def invincibility_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(source=interaction)

        do = DiscordObject(interaction=interaction)

        member_dict = await do.determine_from_target(target=member)
        kwargs = member_dict.get("columns", None)

        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(**kwargs)
            await Invincibility.unrestrict(**kwargs)
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member_dict.get("mention", None)}."
            )
        else:
            Invincibility.remove_invincible_member(**kwargs)
            msg = f"Invincibility has been disabled for {member_dict.get("mention", None)}"

        return await state.end(success=msg)

    # DONE
    @commands.command(name="hero", help="Grant/revoke invincibility.")
    @guild_owner_predicator()
    async def invincibility_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Tag a member or include their ID"
        ),
    ):
        state = StateService(source=ctx)

        do = DiscordObject(ctx=ctx)

        member_dict = await do.determine_from_target(target=member)
        kwargs = member_dict.get("columns", None)

        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(**kwargs)
            await Invincibility.unrestrict(**kwargs)
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member_dict.get("mention", None)}."
            )
        else:
            Invincibility.remove_invincible_member(**kwargs)
            msg = f"Invincibility has been disabled for {member_dict.get("mention", None)}"

        return await state.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerCommands(bot))
