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
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.service.check_service import *
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.member_service import resolve_member
from vyrtuous.service.role_service import resolve_role
from vyrtuous.database.roles.administrator import Administrator, AdministratorRole
from vyrtuous.utils.emojis import get_random_emoji, EMOJIS
from vyrtuous.utils.properties.snowflake import *
from vyrtuous.service.state_service import State
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
        state = State(interaction)
        action = None
        chunk_size = 7
        pages = []
        skipped_members, target_members = [], []
        role_obj = None
        role_snowflakes = []
        try:
            role_obj = await resolve_role(
                ctx_interaction_or_message=interaction, role_str=role
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

        administrators = await Administrator.select_and_role(
            guild_snowflake=interaction.guild.id, role_snowflake=role_obj.id
        )
        administrator_roles = await AdministratorRole.select(
            guild_snowflake=interaction.guild.id
        )
        for administrator_role in administrator_roles:
            await administrator_role.revoke()
        for member in role_obj.members:
            if administrators:
                for administrator in administrators:
                    if role_obj.id not in administrator.role_snowflakes:
                        skipped_members.append(member)
                        continue
                    elif role_obj.id in administrator.role_snowflakes:
                        await administrator.revoke()
                        action = "revoked"
                        target_members.append(member.mention)
        if action != "revoked":
            administrator_role = AdministratorRole(
                guild_snowflake=interaction.guild.id, role_snowflake=role_obj.id
            )
            await administrator_role.grant()
            administrator_roles = await AdministratorRole.select(
                guild_snowflake=interaction.guild.id
            )
            for administrator_role in administrator_roles:
                role_snowflakes.append(administrator_role.role_snowflake)
            for member in role_obj.members:
                administrator = Administrator(
                    guild_snowflake=interaction.guild.id,
                    member_snowflake=member.id,
                    role_snowflakes=role_snowflakes,
                )
                await administrator.grant()
                action = "granted"
                target_members.append(member.mention)

        chunks = [
            target_members[i : i + chunk_size]
            for i in range(0, len(target_members), chunk_size)
        ]
        for index, chunk in enumerate(chunks, start=1):
            embed = discord.Embed(
                title=f"{get_random_emoji()}"
                f"{role_obj.name} Permission Update",
                description=f"Members {action} `Administrator`.",
                color=discord.Color.green(),
            )
            embed.add_field(
                name=f"Members ({len(target_members)})",
                value="\n".join(
                    [
                        m.mention if isinstance(m, discord.Member) else str(m)
                        for m in chunk
                    ]
                ),
                inline=False,
            )
            pages.append(embed)
        if skipped_members:
            chunks = [
                skipped_members[i : i + chunk_size]
                for i in range(0, len(skipped_members), chunk_size)
            ]
            for index, chunk in enumerate(chunks, start=1):
                embed = discord.Embed(
                    title=f"{get_random_emoji()} "
                    f"{role_obj.name} Skipped Members",
                    description=f"Members with {role_obj.mention}",
                    color=discord.Color.red(),
                )
                embed.add_field(
                    name=f"Members ({len(skipped_members)})",
                    value="\n".join(
                        [
                            m.mention if isinstance(m, discord.Member) else str(m)
                            for m in chunk
                        ]
                    ),
                    inline=False,
                )
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f "
                        f"Embed size is too large. Limit the scope."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
        else:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f " f"No members found."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="arole", help="Role -> Administrator.")
    @guild_owner_predicator()
    async def grant_administrator_by_role_text_command(
        self, ctx: commands.Context, role: RoleSnowflake
    ):
        state = State(ctx)
        action = None
        chunk_size = 7
        pages = []
        skipped_members, target_members = [], []
        role_obj = None
        role_snowflakes = []
        try:
            role_obj = await resolve_role(ctx, role)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

        administrators = await Administrator.select_and_role(
            guild_snowflake=ctx.guild.id, role_snowflake=role_obj.id
        )
        administrator_roles = await AdministratorRole.select(
            guild_snowflake=ctx.guild.id
        )
        for administrator_role in administrator_roles:
            await administrator_role.revoke()
        for member in role_obj.members:
            if administrators:
                for administrator in administrators:
                    if role_obj.id not in administrator.role_snowflakes:
                        skipped_members.append(member)
                        continue
                    elif role_obj.id in administrator.role_snowflakes:
                        await administrator.revoke()
                        action = "revoked"
                        target_members.append(member.mention)
        if action != "revoked":
            administrator_role = AdministratorRole(
                guild_snowflake=ctx.guild.id, role_snowflake=role_obj.id
            )
            await administrator_role.grant()
            administrator_roles = await AdministratorRole.select(
                guild_snowflake=ctx.guild.id
            )
            administrator_roles.append(administrator_role)
            for administrator_role in administrator_roles:
                role_snowflakes.append(administrator_role.role_snowflake)
            for member in role_obj.members:
                administrator = Administrator(
                    guild_snowflake=ctx.guild.id,
                    member_snowflake=member.id,
                    role_snowflakes=role_snowflakes,
                )
                await administrator.grant()
                action = "granted"
                target_members.append(member.mention)

        chunks = [
            target_members[i : i + chunk_size]
            for i in range(0, len(target_members), chunk_size)
        ]
        for index, chunk in enumerate(chunks, start=1):
            embed = discord.Embed(
                title=f"{get_random_emoji()} "
                f"{role_obj.name} Permission Update",
                description=f"Members {action} `Administrator`.",
                color=discord.Color.green(),
            )
            embed.add_field(
                name=f"Members ({len(target_members)})",
                value="\n".join(
                    [
                        m.mention if isinstance(m, discord.Member) else str(m)
                        for m in chunk
                    ]
                ),
                inline=False,
            )
            pages.append(embed)
        if skipped_members:
            chunks = [
                skipped_members[i : i + chunk_size]
                for i in range(0, len(skipped_members), chunk_size)
            ]
            for index, chunk in enumerate(chunks, start=1):
                embed = discord.Embed(
                    title=f"{get_random_emoji()} "
                    f"{role_obj.name} Skipped Members",
                    description=f"Members with {role_obj.mention}",
                    color=discord.Color.red(),
                )
                embed.add_field(
                    name=f"Members ({len(skipped_members)})",
                    value="\n".join(
                        [
                            m.mention if isinstance(m, discord.Member) else str(m)
                            for m in chunk
                        ]
                    ),
                    inline=False,
                )
                pages.append(embed)

        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f "
                        f"Embed size is too large. Limit the scope."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
        else:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f " f"No members found."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @app_commands.command(name="hero", description="Grant/revoke invincibility.")
    @app_commands.describe(member="Tag a member or include their ID")
    @guild_owner_predicator()
    async def invincibility_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = State(interaction)
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
                warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
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
        state = State(ctx)
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
                    warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
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
                    warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerCommands(bot))
