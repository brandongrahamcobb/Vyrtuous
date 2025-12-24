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
from vyrtuous.service.discord_message_service import DiscordMessageService
from vyrtuous.service.member_service import MemberService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.vegans import Vegans

class OwnerCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.emoji = Emojis()
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
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
        action = None
        member_obj = await self.member_service.resolve_member(interaction, member)
        if member_obj:
            if member_obj.id == interaction.guild.me.id:
                return
            success = await has_equal_or_higher_role(interaction, member_obj, None)
            if not success:
                return await interaction.response.send_message(content=f"\U0001F6AB You are not allowed to toggle {member_obj.mention}'s role as a developer because they are a higher/or equivalent role than you in {interaction.guild.name}.")
            guilds = []
            developers = await Developer.fetch_by_member(member_snowflake=member_obj.id)
            for developer in developers:
                guilds.append(developer['guild_snowflake'])
            if interaction.guild.id in guilds:
                await Developer.update_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                action = 'revoked'
            else:
                developer = Developer(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                await developer.grant()
                action = 'granted'
        return await interaction.response.send_message(content=f"{self.emoji.get_random_emoji()} Developer access has been {action} in {interaction.guild.name}.")
        
    # DONE
    @commands.command(name='dev', help="Grants/revokes a user's permissions to a bot developer.")
    @is_owner_predicator()
    async def create_developer_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
    ) -> None:
        action = None
        member_obj = await self.member_service.resolve_member(ctx, member)
        if member_obj:
            if member_obj.id == ctx.guild.me.id:
                return
            success = await has_equal_or_higher_role(ctx.message, member_obj, None)
            if not success:
                return await self.handler.send_message(ctx, content=f"\U0001F6AB You are not allowed to toggle {member_obj.mention}'s role as a developer because they are a higher/or equivalent role than you in {ctx.guild.name}.")
            guilds = []
            developers = await Developer.fetch_by_member(member_snowflake=member_obj.id)
            for developer in developers:
                guilds.append(developer['guild_snowflake'])
            if ctx.guild.id in guilds:
                await Developer.update_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                action = 'revoked'
            else:
                developer = Developer(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                await developer.grant()
                action = 'granted'
        await self.handler.send_message(ctx, content=f"{self.emoji.get_random_emoji()} Developer access has been {action} in {ctx.guild.name}.")
        
    # DONE
    @app_commands.command(name='hero', description='Grants/revokes invincibility for a member.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_predicator()
    async def vegan_hero_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake
    ):
        state = None
        member_obj = await self.member_service.resolve_member(interaction, member)
        if member_obj:
            if member_obj.id == interaction.guild.me.id:
                return
            state = Vegans.toggle_state()
            if state:
                Vegans.add_vegan(member_obj.id)
                await Vegans.unrestrict(interaction.guild, member_obj)
            else:
                Vegans.remove_vegan(member_obj.id)
            state = 'ON' if state else 'OFF'
        await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Superhero mode turned {state}.')
           
    # DONE
    @commands.command(name='hero', help='Grants/revokes invincibility for a member.')
    @is_owner_predicator()
    async def vegan_hero_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(description='Tag a member or include their snowflake ID')
    ):
        state = None
        member_obj = await self.member_service.resolve_member(ctx, member)
        if member_obj:
            if member_obj.id == ctx.guild.me.id:
                return
            state = Vegans.toggle_state()
            if state:
                Vegans.add_vegan(member_obj.id)
                await Vegans.unrestrict(ctx.guild, member_obj)
            else:
                Vegans.remove_vegan(member_obj.id)
            state = f'ON' if state else f'OFF'
        await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Superhero mode turned {state}.')

async def setup(bot: DiscordBot):
    await bot.add_cog(OwnerCommands(bot))
