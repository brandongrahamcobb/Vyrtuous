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
    has_equal_or_higher_role,
    not_bot,
    sys_owner_predicator,
)
from vyrtuous.service.logging_service import logger
from vyrtuous.service.messaging.message_service import MessageService
from vyrtuous.service.messaging.state_service import StateService
from vyrtuous.service.resolution.member_service import resolve_member
from vyrtuous.utils.emojis import get_random_emoji


class SystemOwnerCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    # DONE
    @app_commands.command(name="adev", description="Assign developer.")
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
        state = StateService(interaction)
        member_obj = None
        try:
            member_obj = await resolve_member(interaction, member)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                await state.end(error=f"\u274c {str(e).capitalize()}")
        developer = await Developer.select(
            guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id
        )
        if developer:
            developer_log = await DeveloperLog.select(id=reference, resolved=False)
            if developer_log:
                channel = self.bot.get_channel(developer_log.channel_snowflake)
                try:
                    msg = await channel.fetch_message(developer_log.message_snowflake)
                    link = msg.jump_url
                except discord.NotFound:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f Message reference not found: {reference}."
                        )
                    except Exception as e:
                        link = "Unknown message"
                        logger.warning(f"{str(e).capitalize()}")
                        pass
                if developer.member_snowflake in developer_log.developer_snowflakes:
                    await developer_log.unassign(member_snowflake=member_obj.id)
                    try:
                        return await state.end(
                            success=f"{get_random_emoji()} Developer {member_obj.mention} unassigned for issue by {interaction.user.mention}: {link}\n**Notes:** {developer_log.notes}."
                        )
                    except Exception as e:
                        return await state.end(error=f"\u274c {str(e).capitalize()}")
                else:
                    await developer_log.assign(member_snowflake=member_obj.id)
                    try:
                        await member_obj.send(
                            f"{get_random_emoji()} Developer {member_obj.mention} assigned to issue by {interaction.user.mention}: {link}\n**Notes:** {developer_log.notes}"
                        )
                    except Exception as e:
                        logger.warning(f"{str(e).capitalize()}")
                        pass
                    try:
                        return await state.end(
                            success=f"{get_random_emoji()} Developer {member_obj.mention} assigned to issue by {interaction.user.mention}: {link}\n**Notes:** {developer_log.notes}."
                        )
                    except Exception as e:
                        return await state.end(error=f"\u274c {str(e).capitalize()}")
            else:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f Unresolved issue not found for reference: {reference}."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
        else:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f Developer not found for {interaction.guild.name}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="adev", help="Assign developer.")
    @sys_owner_predicator()
    async def toggle_issue_to_developer_text_command(
        self,
        ctx: commands.Context,
        reference: str = commands.parameter(
            default=None, description="Include an issue reference ID"
        ),
        member: MemberSnowflake = commands.parameter(
            default=None, description="Tag a member or include their ID"
        ),
    ):
        state = StateService(ctx)
        member_obj = None
        try:
            member_obj = await resolve_member(ctx, member)
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                await state.end(error=f"\u274c {str(e).capitalize()}")
        developer = await Developer.select(
            guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id
        )
        if developer:
            developer_log = await DeveloperLog.select(id=reference, resolved=False)
            if developer_log:
                channel = self.bot.get_channel(developer_log.channel_snowflake)
                try:
                    msg = await channel.fetch_message(developer_log.message_snowflake)
                    link = msg.jump_url
                except discord.NotFound:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f Message reference not found: {reference}."
                        )
                    except Exception as e:
                        link = "Unknown message"
                        logger.warning(f"{str(e).capitalize()}")
                        pass
                if developer.member_snowflake in developer_log.developer_snowflakes:
                    await developer_log.unassign(member_snowflake=member_obj.id)
                    try:
                        return await state.end(
                            success=f"{get_random_emoji()} Developer {member_obj.mention} unassigned for issue by {ctx.author.mention}: {link}\n**Notes:** {developer_log.notes}."
                        )
                    except Exception as e:
                        return await state.end(error=f"\u274c {str(e).capitalize()}")
                else:
                    await developer_log.assign(member_snowflake=member_obj.id)
                    try:
                        await member_obj.send(
                            f"{get_random_emoji()} Developer {member_obj.mention} assigned for issue by {ctx.author.mention}: {link}\n**Notes:** {developer_log.notes}"
                        )
                    except Exception as e:
                        logger.warning(f"{str(e).capitalize()}")
                        pass
                    try:
                        return await state.end(
                            success=f"{get_random_emoji()} Developer {member_obj.mention} assigned for issue by {ctx.author.mention}: {link}\n**Notes:** {developer_log.notes}."
                        )
                    except Exception as e:
                        return await state.end(error=f"\u274c {str(e).capitalize()}")
            else:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f Unresolved issue not found for reference: {reference}."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
        else:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f Developer not found for {ctx.guild.name}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @app_commands.command(name="dev", description="Grant/revoke devs.")
    @app_commands.describe(member="Tag a member or include their ID")
    @sys_owner_predicator()
    async def create_developer_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(interaction)
        action = None
        member_obj = None
        try:
            member_obj = await resolve_member(interaction, member)
            not_bot(interaction, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(
                interaction,
                channel_snowflake=interaction.channel.id,
                guild_snowflake=interaction.guild.id,
                member_snowflake=member_obj.id,
                sender_snowflake=interaction.user.id,
            )
        except Exception as e:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        developer = await Developer.select(
            guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id
        )
        if developer:
            await Developer.delete(
                guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id
            )
            action = "revoked"
        else:
            developer = Developer(
                guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id
            )
            await developer.create()
            action = "granted"
        try:
            return await state.end(
                success=f"{get_random_emoji()} Developer access for {member_obj.mention} has been {action} in {interaction.guild.name}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="dev", help="Grant/revoke devs.")
    @sys_owner_predicator()
    async def create_developer_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            default=None, description="Tag a member or include their ID"
        ),
    ):
        state = StateService(ctx)
        action = None
        member_obj = None
        try:
            member_obj = await resolve_member(ctx, member)
            not_bot(ctx, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(
                ctx,
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
        developer = await Developer.select(
            guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id
        )
        if developer:
            await Developer.delete(
                guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id
            )
            action = "revoked"
        else:
            developer = Developer(
                guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id
            )
            await developer.create()
            action = "granted"
        try:
            return await state.end(
                success=f"{get_random_emoji()} Developer access for {member_obj.mention} has been {action} in {ctx.guild.name}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")


async def setup(bot: DiscordBot):
    await bot.add_cog(SystemOwnerCommands(bot))
