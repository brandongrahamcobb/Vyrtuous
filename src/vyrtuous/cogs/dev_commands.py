''' dev_commands.py A discord.py cog containing developer commands for the Vyrtuous bot.

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
from typing import Literal, Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.check_service import *
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.role_service import RoleService
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.database import Database
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.state import State
from vyrtuous.utils.temporary_room import TemporaryRoom

class DevCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.handler = MessageService(self.bot, self.bot.db_pool)
        self.member_service = MemberService()
        self.role_service = RoleService()
    
    # DONE
    @app_commands.command(name='backup', description='Creates a backup of the database and uploads it.')
    @is_owner_developer_predicator()
    async def app_backup(
        self,
        interaction: discord.Interaction
    ):
        state = State(interaction)
        db = Database(directory='/app/backups')
        try:
            db.create_backup_directory()
            db.execute_backup()
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @commands.command(name='backup', help='Creates a backup of the database and uploads it.')
    @is_owner_developer_predicator()
    async def backup(
        self,
        ctx: commands.Context
    ):
        state = State(ctx)
        db = Database(directory='/app/backups')
        try:
            db.create_backup_directory()
            db.execute_backup()
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
 
    # DONE
    @app_commands.command(name='clear', description='Removes a specific channel ID from all users, including temp-room associations.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_predicator()
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        await Alias.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        await Coordinator.delete_channel(channel_snowflake=channel_obj.id)
        await Moderator.delete_channel(channel_snowflake=channel_obj.id)
        await TemporaryRoom.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Removed channel `{channel_obj.mention}` from all users\' coordinator and moderator access and deleted all associated records.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    # DONE
    @commands.command(name='clear', help='Removes a specific channel ID from all users, including temp-room associations and all related records.')
    @is_owner_developer_predicator()
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        await Alias.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        await Coordinator.delete_channel(channel_snowflake=channel_obj.id)
        await Moderator.delete_channel(channel_snowflake=channel_obj.id)
        await TemporaryRoom.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Removed channel `{channel_obj.mention}` from all users\' coordinator and moderator access and deleted all associated records.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
        
    @app_commands.command(name='load', description='Loads a cog by name "vyrtuous.cog.<cog_name>.')
    @is_owner_developer_predicator()
    async def load_app_command(self, interaction: discord.Interaction, module: str):
        state = State(interaction)
        try:
            await interaction.response.defer(ephemeral=True)
            await interaction.client.load_extension(module)
        except commands.ExtensionError as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {e}.')
        try:
            return await state.end(f'{self.emoji.get_random_emoji()}  Successfully loaded {module}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
            
    @commands.command(name='load', help='Loads a cog by name "vyrtuous.cog.<cog_name>."')
    @is_owner_developer_predicator()
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.load_extension(module)
        except commands.ExtensionError as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {e}.')
        try:
            return await state.end(f'{self.emoji.get_random_emoji()} Successfully loaded {module}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
    @app_commands.command(name='reload', description='Reloads a cog by name "vyrtuous.cog.<cog_name>".')
    @app_commands.check(at_home)
    @is_owner_developer_predicator()
    async def reload_app_command(self, interaction: discord.Interaction, module: str):
        state = State(interaction)
        try:
            await interaction.response.defer(ephemeral=True)
            await interaction.client.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {e}.')
        try:
            return await state.end(f'{self.emoji.get_random_emoji()} Successfully reloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
            
    @commands.command(name='reload', help='Reloads a cog by name "vyrtuous.cog.<cog_name>".')
    @is_owner_developer_predicator()
    async def reload(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {e}.')
        try:
            return await state.end(f'{self.emoji.get_random_emoji()} Successfully reloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

    @app_commands.command(name='sync', description="Syncs commands to the tree '~', globally '*', clear '^' or general sync.")
    @is_owner_developer_predicator()
    async def sync_app_command(self, interaction: discord.Interaction, spec: Optional[Literal['~', '*', '^']] = None):
        state = State(interaction)
        guilds = interaction.client.guilds
        await interaction.response.defer(ephemeral=True)
        if not guilds:
            if spec == '~':
                synced = await interaction.client.tree.sync(guild=interaction.guild)
            elif spec == '*':
                interaction.client.tree.copy_global_to(guild=interaction.guild)
                synced = await interaction.client.tree.sync(guild=interaction.guild)
            elif spec == '^':
                interaction.client.tree.clear_commands(guild=interaction.guild)
                await interaction.client.tree.sync(guild=interaction.guild)
                synced = []
            else:
                synced = await interaction.client.tree.sync()
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}")
            except Exception as e:   
                return await state.end(error=f'\U0001F3C6 {e}.')
        ret = 0
        for guild in guilds:
            try:
                await interaction.client.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Synced the tree to {ret}/{len(guilds)}.')
        except Exception as e:   
            return await state.end(error=f'\U0001F3C6 {e}.')
        
    @commands.command(name='sync', help="Syncs commands to the tree '~', globally '*', clear '^' or general sync.")
    @is_owner_developer_predicator()
    async def sync_text_command(self, ctx: commands.Context, guilds: commands.Greedy[discord.Object], spec: Optional[Literal['~', '*', '^']] = None):
        state = State(ctx)
        if not guilds:
            if spec == '~':
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '*':
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '^':
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
                synced = []
            else:
                synced = await ctx.bot.tree.sync()
            try:
                return await state.end(success=f"{self.emoji.get_random_emoji()} Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}")
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Synced the tree to {ret}/{len(guilds)}.')
        except Exception as e:   
            return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @app_commands.command(name='trole', description='Marks a role as administrator and syncs all members.')
    @is_owner_developer_predicator()
    async def grant_team_to_role_app_command(
        self,
        interaction: discord.Interaction,
        role: AppRoleSnowflake
    ):
        state = State(interaction)
        try:
            role_obj = await self.role_service.resolve_role(interaction, role)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        for member in role_obj.members:
            administrator = Administrator(guild_snowflake=interaction.guild.id, member_snowflake=member.id, role_snowflake=role_obj.id)
            await administrator.grant()
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Team role granted for members with {role_obj.mention}.')
        except Exception as e:   
            return await state.end(error=f'\U0001F3C6 {e}.')
        
    @commands.command(name='trole', help='Marks a role as administrator and syncs all members.')
    @is_owner_developer_predicator()
    async def grant_team_to_role_text_command(
        self,
        ctx: commands.Context,
        role: RoleSnowflake
    ):
        state = State(ctx)
        try:
            role_obj = await self.role_service.resolve_role(ctx, role)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        for member in role_obj.members:
            administrator = Administrator(guild_snowflake=ctx.guild.id, member_snowflake=member.id, role_snowflake=role_obj.id)
            await administrator.grant()
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Team role granted for members with {role_obj.mention}.')
        except Exception as e:   
            return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @app_commands.command(name='unload', description='Unloads a cog by name "vyrtuous.cog.<cog_name>".')
    @is_owner_developer_predicator()
    async def unload_app_command(self, interaction: discord.Interaction, module: str):
        state = State(interaction)
        try:
            await interaction.response.defer(ephemeral=True)
            await interaction.client.unload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {e}.')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully unloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @commands.command(name='unload', help='Unload a cog by name "vyrtuous.cog.<cog_name>".')
    @is_owner_developer_predicator()
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {e}.')
        try:
            return await state.end(f'{self.emoji.get_random_emoji()} Successfully unloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
                
    @is_owner_developer_predicator()
    async def revoke_administrator_to_role_app_command(
        self,
        interaction: discord.Interaction,
        role: AppRoleSnowflake
    ):
        state = State(interaction)
        try:
            role_obj = await self.role_service.resolve_role(interaction, role)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        for member in role_obj.members:
            administrator = await Administrator.fetch_member(member.id)
            if not administrator:
                continue
            if administrator.role_snowflake == role_obj.id:
                administrator.role_snowflake = None
            if administrator.role_snowflake not in {r.id for r in member.roles}:
                if administrator.guild_snowflake == interaction.guild.id:
                    administrator.guild_snowflake = None
            await Administrator.update_guild_and_role_for_member(
                guild_snowflake=administrator.guild_snowflake,
                member_snowflake=member.id,
                role_snowflake=administrator.role_snowflake
            )
        try:
            return await state.end(f'{self.emoji.get_random_emoji()} Successfully revoked administrator from all members with {role_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

    @commands.command(name='xtrole', help='Revokes a role from administrator and updates all members.')
    @is_owner_developer_predicator()
    async def revoke_administrator_to_role_text_command(
        self,
        ctx: commands.Context,
        role: RoleSnowflake
    ):
        state = State(ctx)
        try:
            role_obj = await self.role_service.resolve_role(ctx, role)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        for member in role_obj.members:
            administrator = await Administrator.fetch_member(member.id)
            if not administrator:
                continue
            if administrator.role_snowflake == role_obj.id:
                administrator.role_snowflake = None
            if administrator.role_snowflake not in {r.id for r in member.roles}:
                if administrator.guild_snowflake == ctx.guild.id:
                    administrator.guild_snowflake = None
            await Administrator.update_guild_and_role_for_member(
                guild_snowflake=administrator.guild_snowflake,
                member_snowflake=member.id,
                role_snowflake=administrator.role_snowflake
            )
        try:
            return await state.end(f'{self.emoji.get_random_emoji()} Successfully revoked administrator from all members with {role_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')    

async def setup(bot: DiscordBot):
    await bot.add_cog(DevCommands(bot))
