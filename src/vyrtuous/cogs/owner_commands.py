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
from vyrtuous.utils.vegans import Vegans

class OwnerCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.emoji = Emojis()
        self.handler = MessageService(self.bot, self.bot.db_pool)
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
            await has_equal_or_higher_role(interaction, channel_snowflake=interaction.channel.id, member_snowflake=member_obj.id)
        except Exception as e:
            await state.end(warning=str(e))
        if member_obj.id == interaction.guild.me.id:
            await state.end(warning="You cannot promote the bot to developer.")
        guilds = []
        developers = await Developer.fetch_by_member(member_snowflake=member_obj.id)
        if developers:
            for developer in developers:
                guilds.append(developer.guild_snowflake)
        if interaction.guild.id in guilds:
            await Developer.update_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            action = 'revoked'
        else:
            developer = Developer(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            await developer.grant()
            action = 'granted'
        try:
            await state.end(success=f"{self.emoji.get_random_emoji()} Developer access for {member_obj.mention} has been {action} in {interaction.guild.name}.")
        except Exception as e:
            await state.end(error=f'\U0001F3C6 {str(e)}')
        
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
            await has_equal_or_higher_role(ctx, channel_snowflake=ctx.channel.id, member_snowflake=member_obj.id)
        except Exception as e:
            await state.end(warning=str(e))
        if member_obj.id == ctx.guild.me.id:
            await state.end(warning="You cannot promote the bot to developer.")
        guilds = []
        developers = await Developer.fetch_by_member(member_snowflake=member_obj.id)
        if developers:
            for developer in developers:
                guilds.append(developer.guild_snowflake)
        if ctx.guild.id in guilds:
            await Developer.update_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            action = 'revoked'
        else:
            developer = Developer(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            await developer.grant()
            action = 'granted'
        try:
            await state.end(success=f"{self.emoji.get_random_emoji()} Developer access for {member_obj.mention} has been {action} in {ctx.guild.name}.")
        except Exception as e:
            await state.end(error=f'\U0001F3C6 {str(e)}')
        
        
    # DONE
    @app_commands.command(name='hero', description='Grants/revokes invincibility for a member.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_predicator()
    async def vegan_hero_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake
    ):
        state = State(interaction)
        enabled = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
        except Exception as e:
            await state.end(warning=str(e))
        if member_obj.id == interaction.guild.me.id:
            await state.end(warning="You cannot promote the bot to a hero.")
        enabled = Vegans.toggle_state()
        if enabled:
            Vegans.add_vegan(member_obj.id)
            await Vegans.unrestrict(interaction.guild, member_obj)
        else:
            Vegans.remove_vegan(member_obj.id)
        enabled = 'ON' if enabled else 'OFF'
        try:
            await self.handler.send_message(interaction, content=f'{self.emoji.get_random_emoji()} Superhero mode turned {enabled} for {member_obj.mention}.')
        except Exception as e:
            await state.end(error=f'\U0001F3C6 {str(e)}')
           
    # DONE
    @commands.command(name='hero', help='Grants/revokes invincibility for a member.')
    @is_owner_predicator()
    async def vegan_hero_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(description='Tag a member or include their snowflake ID')
    ):
        state = State(ctx)
        enabled = None
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
        except Exception as e:
            await state.end(warning=str(e))
        if member_obj.id == ctx.guild.me.id:
            await state.end(warning="You cannot promote the bot to a hero.")
        enabled = Vegans.toggle_enabled()
        if enabled:
            Vegans.add_vegan(member_obj.id)
            await Vegans.unrestrict(ctx.guild, member_obj)
        else:
            Vegans.remove_vegan(member_obj.id)
        enabled = f'ON' if enabled else f'OFF'
        try:
            await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Superhero mode turned {enabled}.')
        except Exception as e:
            await state.end(error=f'\U0001F3C6 {str(e)}')
async def setup(bot: DiscordBot):
    await bot.add_cog(OwnerCommands(bot))
