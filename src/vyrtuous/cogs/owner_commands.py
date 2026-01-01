''' owner_commands.py A discord.py cog containing owner-only commands for the Vyrtuous bot.

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
'''
from discord import app_commands
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.service.check_service import *
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.member_service import MemberService
from vyrtuous.utils.developer import Developer
from vyrtuous.utils.developer_log import DeveloperLog
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.state import State
from vyrtuous.utils.invincibility import Invincibility

class OwnerCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.emoji = Emojis()
        self.message_service = MessageService(self.bot, self.bot.db_pool)
        self.member_service = MemberService()

    @app_commands.command(name='dish', description="Assigns or unassigns a developer to a pending issue.")
    @app_commands.describe(
        reference='Include an issue reference ID',
        member='Tag a member or include their snowflake ID',
    )
    @is_owner_predicator()
    async def toggle_issue_to_developer_app_command(
        self,
        interaction: discord.Interaction,
        reference: Optional[str],
        member: AppMemberSnowflake,
    ):
        state = State(interaction)
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                await state.end(error=f'\u274C {str(e).capitalize()}')
        developers = await Developer.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
        if developers:
            developer_log = await DeveloperLog.fetch_unresolved_by_reference(reference)
            for developer_snowflake in developer_log.developer_snowflakes:
                if member_obj.id == developer_snowflake:
                    await developer_log.unassign(member_snowflake=member_obj.id)
                    try:
                        return await state.end(success=f'{self.emoji.get_random_emoji()} Developer unassigned for issue: {msg.link_url}\n**Notes:** {developer_log.notes}.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
            if developer_log:
                developer_log.assign(member_snowflake=member_obj.id)
                channel = self.bot.get_channel(developer_log.channel_snowflake)
                try:
                    msg = await channel.fetch_message(developer_log.message_snowflake)
                    await member_obj.send('{self.emoji.get_random_emoji()} Developer assigned to issue: {msg.link_url}\n**Notes:** {developer_log.notes}')
                    return await state.end(success=f'{self.emoji.get_random_emoji()} Developer assigned to issue: {msg.link_url}\n**Notes:** {developer_log.notes}.')
                except Exception as e:
                    try:
                        return await state.end(success=f'{self.emoji.get_random_emoji()} Developer assigned to issue: {msg.link_url}\n**Notes:** {developer_log.notes}.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
            else:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Unresolved issue not found for reference: {reference}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Developer not found for {interaction.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    @commands.command(name='dish', help="Assigns or unassigns a developer to a pending issue.")
    @is_owner_predicator()
    async def toggle_issue_to_developer_text_command(
        self,
        ctx: commands.Context,
        reference: Optional[str] = commands.parameter(default=None, description='Include an issue reference ID'),
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
    ):
        state = State(ctx)
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                await state.end(error=f'\u274C {str(e).capitalize()}')
        developers = await Developer.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
        if developers:
            developer_log = await DeveloperLog.fetch_unresolved_by_reference(reference)
            for developer_snowflake in developer_log.developer_snowflakes:
                if member_obj.id == developer_snowflake:
                    await developer_log.unassign(member_snowflake=member_obj.id)
                    try:
                        return await state.end(success=f'{self.emoji.get_random_emoji()} Developer unassigned for issue: {msg.link_url}\n**Notes:** {developer_log.notes}.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
            if developer_log:
                developer_log.assign(member_snowflake=member_obj.id)
                channel = self.bot.get_channel(developer_log.channel_snowflake)
                try:
                    msg = await channel.fetch_message(developer_log.message_snowflake)
                    await member_obj.send('{self.emoji.get_random_emoji()} Developer assigned to issue: {msg.link_url}\n**Notes:** {developer_log.notes}')
                    return await state.end(success=f'{self.emoji.get_random_emoji()} Developer assigned to issue: {msg.link_url}\n**Notes:** {developer_log.notes}.')
                except Exception as e:
                    try:
                        return await state.end(success=f'{self.emoji.get_random_emoji()} Developer assigned to issue: {msg.link_url}\n**Notes:** {developer_log.notes}.')
                    except Exception as e:
                        return await state.end(error=f'\u274C {str(e).capitalize()}')
            else:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F Unresolved issue not found for reference: {reference}.')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Developer not found for {ctx.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='dev', description="Grants/revokes a user's permissions to a bot developer.")
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_predicator()
    async def create_developer_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake
    ):
        state = State(interaction)
        action = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
            check_not_self(interaction, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(interaction, channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, sender_snowflake=interaction.user.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        developer = await Developer.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
        if developer:
            await developer.revoke()
            action = 'revoked'
        else:
            developer = Developer(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            await developer.grant()
            action = 'granted'
        try:
            return await state.end(success=f"{self.emoji.get_random_emoji()} Developer access for {member_obj.mention} has been {action} in {interaction.guild.name}.")
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='dev', help="Grants/revokes a user's permissions to a bot developer.")
    @is_owner_predicator()
    async def create_developer_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
    ):
        state = State(ctx)
        action = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
            check_not_self(ctx, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(ctx, channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, sender_snowflake=ctx.author.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        developer = await Developer.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
        if developer:
            await developer.revoke()
            action = 'revoked'
        else:
            developer = Developer(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            await developer.grant()
            action = 'granted'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Developer access for {member_obj.mention} has been {action} in {ctx.guild.name}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
        
    # DONE
    @app_commands.command(name='hero', description='Grants/revokes invincibility for a member.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_predicator()
    async def invincibility_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake
    ):
        state = State(interaction)
        enabled = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
            check_not_self(interaction, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(interaction, channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, sender_snowflake=interaction.user.id)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(member_snowflake=member_obj.id)
            await Invincibility.unrestrict(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            msg = f'All moderation events have been forgiven and invincibility has been enabled for {member_obj.mention}.'
        else:
            Invincibility.remove_invincible_member(member_snowflake=member_obj.id)
            msg = f'Invincibiility has been disabled for {member_obj.mention}'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
           
    # DONE
    @commands.command(name='hero', help='Grants/revokes invincibility for a member.')
    @is_owner_predicator()
    async def invincibility_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(description='Tag a member or include their snowflake ID')
    ):
        state = State(ctx)
        enabled = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
            check_not_self(ctx, member_snowflake=member_obj.id)
            await has_equal_or_higher_role(ctx, channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, sender_snowflake=ctx.author.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(member_snowflake=member_obj.id)
            await Invincibility.unrestrict(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            msg = f'All moderation events have been forgiven and invincibility has been enabled for {member_obj.mention}.'
        else:
            Invincibility.remove_invincible_member(member_snowflake=member_obj.id)
            msg = f'Invincibiility has been disabled for {member_obj.mention}'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
async def setup(bot: DiscordBot):
    await bot.add_cog(OwnerCommands(bot))
