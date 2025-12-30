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
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.service.check_service import *
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.member_service import MemberService
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
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}')
            except Exception as e:
                await state.end(error=f'\U0001F3C6 {e}')
        guild_snowflakes = await Developer.fetch_guilds_by_member(member_snowflake=member_obj.id)
        if interaction.guild.id in guild_snowflakes:
            await Developer.delete_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            action = 'revoked'
        else:
            developer = Developer(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            await developer.grant()
            action = 'granted'
        try:
            await state.end(success=f"{self.emoji.get_random_emoji()} Developer access for {member_obj.mention} has been {action} in {interaction.guild.name}.")
        except Exception as e:
            await state.end(error=f'\U0001F3C6 {e}')
        
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
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}')
            except Exception as e:
                await state.end(error=f'\U0001F3C6 {e}')
        guild_snowflakes = await Developer.fetch_guilds_by_member(member_snowflake=member_obj.id)
        if ctx.guild.id in guild_snowflakes:
            await Developer.delete_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            action = 'revoked'
        else:
            developer = Developer(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            await developer.grant()
            action = 'granted'
        try:
            await state.end(success=f'{self.emoji.get_random_emoji()} Developer access for {member_obj.mention} has been {action} in {ctx.guild.name}.')
        except Exception as e:
            await state.end(error=f'\U0001F3C6 {e}')
        
        
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
            await has_equal_or_higher_role(interaction, channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, sender_snowflake=interaction.author.id)
        except Exception as e:
            await state.end(warning=f'\U000026A0\U0000FE0F {e}')
        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(member_obj.id)
            await Invincibility.unrestrict(interaction.guild, member_obj)
            msg = f'All moderation events have been forgiven and invincibility has been enabled for {member_obj.mention}.'
        else:
            Invincibility.remove_invincible_member(member_obj.id)
            msg = f'Invincibiility has been disabled for {member_obj.mention}'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
           
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
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(member_obj.id)
            await Invincibility.unrestrict(ctx.guild, member_obj)
            msg = f'All moderation events have been forgiven and invincibility has been enabled for {member_obj.mention}.'
        else:
            Invincibility.remove_invincible_member(member_obj.id)
            msg = f'Invincibiility has been disabled for {member_obj.mention}'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}')
        
async def setup(bot: DiscordBot):
    await bot.add_cog(OwnerCommands(bot))
