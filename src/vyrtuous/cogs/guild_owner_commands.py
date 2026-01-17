"""guild_owner_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

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
from vyrtuous.database.roles.administrator import Administrator, AdministratorRole
from vyrtuous.properties.snowflake import (
    AppMemberSnowflake,
    AppRoleSnowflake,
    MemberSnowflake,
    RoleSnowflake,
)
from vyrtuous.service.check_service import (
    guild_owner_predicator,
    has_equal_or_higher_role,
    not_bot,
)
from vyrtuous.service.messaging.message_service import MessageService
from vyrtuous.service.messaging.state_service import StateService
from vyrtuous.service.resolution.member_service import resolve_member
from vyrtuous.service.resolution.role_service import resolve_role
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.invincibility import Invincibility


class GuildOwnerCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    # DONE
    @app_commands.command(name="arole", description="Role -> Administrator.")
    @guild_owner_predicator()
    async def grant_administrator_by_role_app_command(
        self, interaction: discord.Interaction, role: AppRoleSnowflake
    ):
        state = StateService(interaction)

        try:
            role_obj = await resolve_role(ctx_interaction_or_message=interaction, role_str=role)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

        administrators = await Administrator.select(role_snowflakes=[role_obj.id])
        administrator_roles = await AdministratorRole.select(role_snowflake=role_obj.id)
        title = f"{get_random_emoji()} Administrators and Roles"

        chunk_size, field_count, pages = 7, 0, []
        action = None

        if administrator_roles:
            for administrator_role in administrator_roles:
                await AdministratorRole.delete(
                    guild_snowflake=administrator_role.guild_snowflake,
                    role_snowflake=administrator_role.role_snowflake,
                )
            action = "revoked"
        else:
            action = "granted"
        for administrator in administrators:
            revoked_members = {}
            member = interaction.guild.get_member(administrator.member_snowflake)
            await Administrator.delete(
                guild_snowflake=administrator.guild_snowflake,
                member_snowflake=administrator.member_snowflake,
            )
            revoked_members.setdefault(administrator.guild_snowflake, {}).setdefault(
                role_obj.id, []
            ).append(member)
        if action != "revoked":
            granted_members = {}
            granted_members.setdefault(role_obj.guild.id, {})[role_obj.id] = []
            administrator_role = AdministratorRole(
                guild_snowflake=role_obj.guild.id,
                role_snowflake=role_obj.id,
            )
            await administrator_role.create()
            role_snowflakes = [role_obj.id]
            for member in role_obj.members:
                administrator = Administrator(
                    guild_snowflake=role_obj.guild.id,
                    member_snowflake=member.id,
                    role_snowflakes=role_snowflakes,
                )
                await administrator.create()
            granted_members[role_obj.guild.id][role_obj.id].append(member)

        embed = discord.Embed(title=title, description=f'`{role_obj.name}` was `{action}`.', color=discord.Color.red() if action == 'revoked' else discord.Color.green())
        embed.add_field(name='Role ID', value=str(role_obj.id), inline=False)
        embed.add_field(name='Guild', value=role_obj.guild.name, inline=False)
        pages.append(embed)
        members = revoked_members.get(role_obj.guild.id, {}).get(role_obj.id, []) if action == 'revoked' else granted_members.get(role_obj.guild.id, {}).get(role_obj.id, [])
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
            embed = discord.Embed(title=f'Members {action.capitalize()}', color=discord.Color.red() if action == 'revoked' else discord.Color.green())
            for member in chunk:
                embed.add_field(name=f'{field_count}. {member}', value=f'{member.mention} ({member.id})', inline=False)
                field_count += 1
            embed.set_footer(text=f'Page {page_number}')
            pages.append(embed)
            page_number += 1

        await StateService.send_pages(obj=AdministratorRole, pages=pages, state=state)

    # DONE
    @commands.command(name="arole", help="Role -> Administrator.")
    @guild_owner_predicator()
    async def grant_administrator_by_role_text_command(
        self, ctx: commands.Context, role: RoleSnowflake
    ):
        state = StateService(ctx)

        try:
            role_obj = await resolve_role(ctx_interaction_or_message=ctx, role_str=role)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

        administrators = await Administrator.select(role_snowflakes=[role_obj.id])
        administrator_roles = await AdministratorRole.select(role_snowflake=role_obj.id)
        title = f"{get_random_emoji()} Administrators and Roles"

        chunk_size, field_count, pages = 7, 0, []
        action = None

        if administrator_roles:
            for administrator_role in administrator_roles:
                await AdministratorRole.delete(
                    guild_snowflake=administrator_role.guild_snowflake,
                    role_snowflake=administrator_role.role_snowflake,
                )
            action = "revoked"
        else:
            action = "granted"
        for administrator in administrators:
            revoked_members = {}
            member = ctx.guild.get_member(administrator.member_snowflake)
            await Administrator.delete(
                guild_snowflake=administrator.guild_snowflake,
                member_snowflake=administrator.member_snowflake,
            )
            revoked_members.setdefault(administrator.guild_snowflake, {}).setdefault(
                role_obj.id, []
            ).append(member)
        if action != "revoked":
            granted_members = {}
            granted_members.setdefault(role_obj.guild.id, {})[role_obj.id] = []
            administrator_role = AdministratorRole(
                guild_snowflake=role_obj.guild.id,
                role_snowflake=role_obj.id,
            )
            await administrator_role.create()
            role_snowflakes = [role_obj.id]
            for member in role_obj.members:
                administrator = Administrator(
                    guild_snowflake=role_obj.guild.id,
                    member_snowflake=member.id,
                    role_snowflakes=role_snowflakes,
                )
                await administrator.create()
            granted_members[role_obj.guild.id][role_obj.id].append(member)

        embed = discord.Embed(title=title, description=f'`{role_obj.name}` was `{action}`.', color=discord.Color.red() if action == 'revoked' else discord.Color.green())
        embed.add_field(name='Role ID', value=str(role_obj.id), inline=False)
        embed.add_field(name='Guild', value=role_obj.guild.name, inline=False)
        pages.append(embed)
        members = revoked_members.get(role_obj.guild.id, {}).get(role_obj.id, []) if action == 'revoked' else granted_members.get(role_obj.guild.id, {}).get(role_obj.id, [])
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
            embed = discord.Embed(title=f'Members {action.capitalize()}', color=discord.Color.red() if action == 'revoked' else discord.Color.green())
            for member in chunk:
                embed.add_field(name=f'{field_count}. {member}', value=f'{member.mention} ({member.id})', inline=False)
                field_count += 1
            embed.set_footer(text=f'Page {page_number}')
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
        state = StateService(interaction)
        enabled = None
        member_obj = None
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=interaction, member_str=member
            )
            not_bot(
                ctx_interaction_or_message=interaction, member_snowflake=member_obj.id
            )
            await has_equal_or_higher_role(
                ctx_interaction_or_message=interaction,
                channel_snowflake=interaction.channel.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=interaction.user.id,
            )
        except Exception as e:
            return await state.end(
                warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
            )
        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(member_snowflake=member_obj.id)
            await Invincibility.unrestrict(
                guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id
            )
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member_obj.mention}."
            )
        else:
            Invincibility.remove_invincible_member(member_snowflake=member_obj.id)
            msg = f"Invincibility has been disabled for {member_obj.mention}"
        try:
            return await state.end(success=f"{get_random_emoji()} {msg}")
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

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
        state = StateService(ctx)
        enabled = None
        member_obj = None
        try:
            member_obj = await resolve_member(
                ctx_interaction_or_message=ctx, member_str=member
            )
            not_bot(ctx_interaction_or_message=ctx, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(
                ctx_interaction_or_message=ctx,
                channel_snowflake=ctx.channel.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=ctx.author.id,
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(member_snowflake=member_obj.id)
            await Invincibility.unrestrict(
                guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id
            )
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member_obj.mention}."
            )
        else:
            Invincibility.remove_invincible_member(member_snowflake=member_obj.id)
            msg = f"Invincibility has been disabled for {member_obj.mention}"
        try:
            return await state.end(success=f"{get_random_emoji()} {msg}")
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerCommands(bot))
