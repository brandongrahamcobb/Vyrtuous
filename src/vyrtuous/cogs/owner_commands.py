''' dev_commands.py

    Copyright (C) 2024  github.com/brandongrahamcobb

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
from typing import Literal, Optional

from vyrtuous.inc.helpers import *

from vyrtuous.service.check_service import *
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.discord_message_service import DiscordMessageService
from vyrtuous.utils.vegans import Vegans

class OwnerCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.member_service = MemberService()

    # DONE
    @app_commands.command(name='dev', description='Elevates a user\'s permissions to a bot developer.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_predicator()
    async def create_developer_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(interaction, member)
        if member_obj.bot:
            return await interaction.response.send_message(content='\U0001F6AB You cannot make the bot a developer.')
        if not member_obj:
            return await interaction.response.send_message(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        success = await has_equal_or_higher_role(interaction.message, member_obj)
        if not success:
            return await interaction.response.send_message(content=f"\U0001F6AB You are not allowed to toggle {member_obj.mention}'s role as a developer because they are a higher/or equivalent role than you in {channel_obj.mention}.", allowed_mentions=discord.AllowedMentions.none())
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT developer_guild_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            action = None
            if row and interaction.guild.id in row['developer_guild_ids']:
                await conn.execute('''
                    UPDATE users
                    SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                        updated_at = NOW()
                    WHERE discord_snowflake = $1
                ''', member_obj.id, interaction.guild.id)
                action = 'revoked'
            else:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, developer_guild_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET developer_guild_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(u.developer_guild_ids || EXCLUDED.developer_guild_ids)
                        )
                        FROM users u WHERE u.discord_snowflake = EXCLUDED.discord_snowflake
                    ),
                    updated_at = NOW()
                ''', member_obj.id, interaction.guild.id)
                action = 'granted'
        return await interaction.response.send_message(content=f"{self.emoji.get_random_emoji()} {member_obj.mention}'s developer accessed has been {action} in {interaction.guild.name}.", allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @commands.command(name='dev', help='Elevates a user\'s permissions to a bot developer.')
    @is_owner_predicator()
    async def create_developer_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a developer.')
        success = await has_equal_or_higher_role(ctx.message, member_obj)
        if not success:
            return await self.handler.send_message(ctx, content=f"\U0001F6AB You are not allowed to toggle {member_obj.mention}'s role as a developer because they are a higher/or equivalent role than you in {channel_obj.mention}.", allowed_mentions=discord.AllowedMentions.none())
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT developer_guild_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            action = None
            if row and interaction.guild.id in row['developer_guild_ids']:
                await conn.execute('''
                    UPDATE users
                    SET developer_guild_ids = array_remove(developer_guild_ids, $2),
                        updated_at = NOW()
                    WHERE discord_snowflake = $1
                ''', member_obj.id, ctx.guild.id)
                action = 'revoked'
            else:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, developer_guild_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET developer_guild_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(u.developer_guild_ids || EXCLUDED.developer_guild_ids)
                        )
                        FROM users u WHERE u.discord_snowflake = EXCLUDED.discord_snowflake
                    ),
                    updated_at = NOW()
                ''', member_obj.id, ctx.guild.id)
                action = 'granted'
        await self.handler.send_message(ctx, content=f"{self.emoji.get_random_emoji()} {member_obj.mention}'s developer accessed has been {action} in {ctx.guild.name}.", allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @app_commands.command(name='hero', description='Toggle the special feature for a member.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_predicator()
    async def vegan_hero_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        member_obj = await self.member_service.resolve_member(interaction, member)
        if not member_obj:
            return await interaction.response.send_message(content='\U0001F6AB Could not resolve the member.')
        if member_obj.bot:
            return await interaction.response.send_message(content='\U0001F6AB You cannot make the bot a superhero.')
        vegans = Vegans.get_vegans()
        vegans.state = not vegans.state
        if vegans.state:
            vegans.add(member_obj.id)
        else:
            vegans.discard(member_obj.id)
            for channel in interaction.guild.channels:
                await self.unrestrict(interaction.guild, member_obj)
        state = 'ON' if vegans.state else 'OFF'
        await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Superhero mode turned {state}')
           
    # DONE
    @commands.command(name='hero', hidden=True)
    @is_owner_predicator()
    async def vegan_hero_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(description='Tag a member or include their snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content='\U0001F6AB Could not resolve the member.')
        if member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a superhero.')
        vegans = Vegans.get_vegans()
        vegans.state = not vegans.state
        if vegans.state:
            vegans.add(member_obj.id)
        else:
            vegans.discard(member_obj.id)
            for channel in ctx.guild.channels:
                await self.unrestrict(ctx.guild, member_obj)
        state = f'ON' if vegans.state else f'OFF'
        await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Superhero modeturned {state}.')

async def setup(bot: DiscordBot):
    await bot.add_cog(OwnerCommands(bot))
