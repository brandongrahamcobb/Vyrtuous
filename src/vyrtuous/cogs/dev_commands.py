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
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.database import Database
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.flag import Flag
from vyrtuous.utils.invincibility import Invincibility
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.state import State
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.vegan import Vegan
from vyrtuous.utils.voice_mute import VoiceMute

class DevCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.message_service = MessageService(self.bot, self.bot.db_pool)
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
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    # DONE
    @commands.command(name='backup', help='Creates a backup of the database and uploads it.')
    @is_owner_developer_predicator()
    async def text_backup(
        self,
        ctx: commands.Context
    ):
        state = State(ctx)
        db = Database(directory='/app/backups')
        try:
            db.create_backup_directory()
            db.execute_backup()
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
 
    # DONE
    @app_commands.command(name='clear', description='Removes a specific channel ID from all users, including temp-room associations.')
    @app_commands.describe(scope='Tag a channel/member or include the snowflake ID')
    @is_owner_developer_predicator()
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        scope: str
    ):
        state = State(interaction)
        channel_obj = None
        highest_role = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, scope)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        except Exception as e:
            try:
                member_obj = await self.member_service.resolve_member(interaction, scope)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        if channel_obj and highest_role in ('Owner', 'Developer'):
            await Alias.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await Ban.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await Coordinator.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await Flag.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await Moderator.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await TemporaryRoom.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await TextMute.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await Vegan.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await VoiceMute.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            msg = f'Deleted all associated moderation actions and roles for {channel_obj.mention}.'
        elif member_obj:
            await Invincibility.unrestrict(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            msg = f'Deleted all associated moderation actions on {member_obj.mention}.'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='clear', help='Removes a specific channel ID from all users, including temp-room associations and all related records.')
    @is_owner_developer_predicator()
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        scope: str = commands.parameter(default=None, description='Tag a channel, a member or include its the snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        highest_role = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, scope)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        except Exception as e:
            try:
                member_obj = await self.member_service.resolve_member(ctx, scope)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
        if channel_obj and highest_role in ('Owner', 'Developer'):
            await Alias.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await Ban.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await Coordinator.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await Flag.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await Moderator.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await TemporaryRoom.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await TextMute.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await Vegan.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await VoiceMute.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            msg = f'Deleted all associated moderation actions and roles for {channel_obj.mention}.'
        elif member_obj:
            await Invincibility.unrestrict(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            msg = f'Deleted all associated moderation actions on {member_obj.mention}.'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    @app_commands.command(name='cogs', description='List cogs.')
    @is_owner_developer_predicator()
    async def list_cogs_app_command(self, interaction: discord.Interaction):
        state = State(interaction)
        loaded, not_loaded = [], []
        embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} Cogs for {interaction.guild.me.name}', color=discord.Color.blurple())
        for cog in sorted(DISCORD_COGS_CLASSES):
            if cog in self.bot.cogs:
                loaded.append(cog)
            else:
                not_loaded.append(cog)
        if loaded:
            embed.add_field(name='Loaded', value='\n'.join(loaded), inline=False)
        if not_loaded:
            embed.add_field(name='Not Loaded', value='\n'.join(not_loaded), inline=False)
        if not loaded and not not_loaded:
            embed.add_field(name='No cogs available.', inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    @commands.command(name='cogs', help='List cogs.')
    @is_owner_developer_predicator()
    async def list_cogs_text_command(self, ctx: commands.Context):
        state = State(ctx)
        loaded, not_loaded = [], []
        embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} Cogs for {ctx.guild.me.name}', color=discord.Color.blurple())
        for cog in sorted(DISCORD_COGS_CLASSES):
            if cog in self.bot.cogs:
                loaded.append(cog)
            else:
                not_loaded.append(cog)
        if loaded:
            embed.add_field(name='Loaded', value='\n'.join(loaded), inline=False)
        if not_loaded:
            embed.add_field(name='Not Loaded', value='\n'.join(not_loaded), inline=False)
        if not loaded and not not_loaded:
            embed.add_field(name='No cogs available.', inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
                
    # DONE
    @app_commands.command(name='load', description="Loads a cog by name 'vyrtuous.cog.<cog_name>.'")
    @is_owner_developer_predicator()
    async def load_app_command(self, interaction: discord.Interaction, module: str):
        state = State(interaction)
        try:
            await interaction.response.defer(ephemeral=True)
            await interaction.client.load_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully loaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE    
    @commands.command(name='load', help="Loads a cog by name 'vyrtuous.cog.<cog_name>.'")
    @is_owner_developer_predicator()
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.load_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully loaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
    
    # DONE
    @app_commands.command(name='ping', description='Ping the bot!')
    @is_owner_developer_predicator()
    async def ping_app_command(
        self,
        interaction: discord.Interaction
    ):
        state = State(interaction)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Pong!')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='ping', description='Ping the bot!')
    @is_owner_developer_predicator()
    async def ping_text_command(
        self,
        ctx: commands.Context
    ):
        state = State(ctx)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Pong!')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')   
    # DONE
    @app_commands.command(name='reload', description="Reloads a cog by name 'vyrtuous.cog.<cog_name>'.")
    @app_commands.check(at_home)
    @is_owner_developer_predicator()
    async def reload_app_command(self, interaction: discord.Interaction, module: str):
        state = State(interaction)
        try:
            await interaction.response.defer(ephemeral=True)
            await interaction.client.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully reloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
            
    # DONE
    @commands.command(name='reload', help="Reloads a cog by name 'vyrtuous.cog.<cog_name>'.")
    @is_owner_developer_predicator()
    async def reload_text_command(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully reloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='sync', description="Syncs commands to the tree '~', globally '*', clear '^' or general sync.")
    @is_owner_developer_predicator()
    async def sync_app_command(self, interaction: discord.Interaction, spec: Optional[Literal['~', '*', '^']] = None):
        state = State(interaction)
        guilds = interaction.client.guilds
        synced = []
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
            else:
                synced = await interaction.client.tree.sync()
            try:
                if spec is None:
                    msg = f'Synced {len(synced)} commands globally.'
                else:
                    msg = f'Synced {len(synced)} commands to the current server.'
                return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
            except Exception as e:   
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
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
            return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='sync', help="Syncs commands to the tree '~', globally '*', clear '^' or general sync.")
    @is_owner_developer_predicator()
    async def sync_text_command(self, ctx: commands.Context, guilds: commands.Greedy[discord.Object], spec: Optional[Literal['~', '*', '^']] = None):
        state = State(ctx)
        synced = []
        if not guilds:
            if spec == '~':
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '*':
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '^':
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
            else:
                synced = await ctx.bot.tree.sync()
            try:
                if spec is None:
                    msg = f'Synced {len(synced)} commands globally.'
                else:
                    msg = f'Synced {len(synced)} commands to the current server.'
                return await state.end(success=f'{self.emoji.get_random_emoji()} {msg}')
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
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
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @app_commands.command(name='trole', description='Marks a role as administrator and syncs all members.')
    @is_owner_developer_predicator()
    async def grant_administrator_by_role_app_command(
        self,
        interaction: discord.Interaction,
        role: AppRoleSnowflake
    ):
        state = State(interaction)
        pages = []
        skipped_members, target_members = [], []
        role_obj = None
        try:
            role_obj = await self.role_service.resolve_role(interaction, role)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e: 
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        administrators = Administrator.fetch_by_guild(guild_snowflake=interaction.guild.id)
        for member in role_obj.members:
            for administrator in administrators:
                if role_obj.id in administrator.snowflakes:
                    skipped_members(member)
                    continue
                else:
                    administrator = Administrator(guild_snowflake=interaction.guild.id, member_snowflake=member.id, role_snowflake=role_obj.id)
                    await administrator.grant()
                    target_members.append(member.mention)
        chunks = [target_members[i:i + 18] for i in range(0, len(target_members), 18)]
        for index, chunk in enumerate(chunks, start=1):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {role_obj.name} Permission Update',
                description=f'Members granted `Administrator`.',
                color=discord.Color.green()
            )
            embed.add_field(name=f'Member(s) ({len(target_members)})', value='\n'.join(chunk), inline=False)
            pages.append(embed)
        if skipped_members:
            chunks = [skipped_members[i:i + 18] for i in range(0, len(skipped_members), 18)]
            for index, chunk in enumerate(chunks, start=1):
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} {role_obj.name} Skipped Members',
                    description=f'Members with {role_obj.mention}',
                    color=discord.Color.red()
                )
                embed.add_field(name=f'Member(s) ({len(skipped_members)})', value='\n'.join(chunk), inline=False)
                pages.append(embed)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No members found with the role {role_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        
    # DONE
    @commands.command(name='trole', help='Marks a role as administrator and syncs all members.')
    @is_owner_developer_predicator()
    async def grant_administrator_by_role_text_command(
        self,
        ctx: commands.Context,
        role: RoleSnowflake
    ):
        state = State(ctx)
        pages = []
        skipped_members, target_members = [], []
        role_obj = None
        try:
            role_obj = await self.role_service.resolve_role(ctx, role)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e: 
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        administrators = Administrator.fetch_by_guild(guild_snowflake=ctx.guild.id)
        for member in role_obj.members:
            for administrator in administrators:
                if role_obj.id in administrator.snowflakes:
                    skipped_members(member)
                    continue
                else:
                    administrator = Administrator(guild_snowflake=ctx.guild.id, member_snowflake=member.id, role_snowflake=role_obj.id)
                    await administrator.grant()
                    target_members.append(member.mention)
        chunks = [target_members[i:i + 18] for i in range(0, len(target_members), 18)]
        for index, chunk in enumerate(chunks, start=1):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {role_obj.name} Permission Update',
                description=f'Members granted `Administrator`.',
                color=discord.Color.green()
            )
            embed.add_field(name=f'Member(s) ({len(target_members)})', value='\n'.join(chunk), inline=False)
            pages.append(embed)
        if skipped_members:
            chunks = [skipped_members[i:i + 18] for i in range(0, len(skipped_members), 18)]
            for index, chunk in enumerate(chunks, start=1):
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} {role_obj.name} Skipped Members',
                    description=f'Members with {role_obj.mention}',
                    color=discord.Color.red()
                )
                embed.add_field(name=f'Member(s) ({len(skipped_members)})', value='\n'.join(chunk), inline=False)
                pages.append(embed)

    # DONE
    @app_commands.command(name='unload', description="Unloads a cog by name 'vyrtuous.cog.<cog_name>'.")
    @is_owner_developer_predicator()
    async def unload_app_command(self, interaction: discord.Interaction, module: str):
        state = State(interaction)
        try:
            await interaction.response.defer(ephemeral=True)
            await interaction.client.unload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully unloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

    # DONE
    @commands.command(name='unload', help="Unload a cog by name 'vyrtuous.cog.<cog_name>'.")
    @is_owner_developer_predicator()
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        state = State(ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e.__class__.__name__}: {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully unloaded {module}.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')
                
    # DONE
    @app_commands.command(name='xtrole', description='Revokes a role from administrator and updates all members.')
    @is_owner_developer_predicator()
    async def revoke_administrator_by_role_app_command(
        self,
        interaction: discord.Interaction,
        role: AppRoleSnowflake
    ):
        state = State(interaction)
        chunk_size = 18
        members_revoked, pages = [], []
        role_obj = None
        try:
            role_obj = await self.role_service.resolve_role(interaction, role)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        for member in role_obj.members:
            administrator = await Administrator.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member.id)
            if not administrator:
                continue
            if role_obj.id in administrator.role_snowflakes:
                await administrator.update_by_removed_role(role_snowflake=role_obj.id)
                members_revoked.append(member.mention)
        member_sets = [members_revoked[i:i + chunk_size] for i in range(0, len(members_revoked), chunk_size)]
        for index, page in enumerate(member_sets, start=1):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {role_obj.name} Permissions Update',
                description=f'Members revoked `Administrator`',
                color=discord.Color.red()
            )
            embed.add_field(
                name=f'Member(s) ({len(page)})',
                value='\n'.join(page),
                inline=False
            )
            pages.append(embed)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No admins found in {interaction.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')    

    # DONE
    @commands.command(name='xtrole', help='Revokes a role from administrator and updates all members.')
    @is_owner_developer_predicator()
    async def revoke_administrator_by_role_text_command(
        self,
        ctx: commands.Context,
        role: RoleSnowflake
    ):
        state = State(ctx)
        chunk_size = 18
        members_revoked, pages = [], []
        role_obj = None
        try:
            role_obj = await self.role_service.resolve_role(ctx, role)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e).capitalize()}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        for member in role_obj.members:
            administrator = await Administrator.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member.id)
            if not administrator:
                continue
            if role_obj.id in administrator.role_snowflakes:
                await administrator.update_by_removed_role(role_snowflake=role_obj.id)
                members_revoked.append(member.mention)
        member_sets = [members_revoked[i:i + chunk_size] for i in range(0, len(members_revoked), chunk_size)]
        for index, page in enumerate(member_sets, start=1):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} {role_obj.name} Permissions Update',
                description=f'Members revoked `Administrator`',
                color=discord.Color.red()
            )
            embed.add_field(
                name=f'Member(s) ({len(page)})',
                value='\n'.join(page),
                inline=False
            )
            pages.append(embed)
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        else:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No admins found in {ctx.guild.name}.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')    

async def setup(bot: DiscordBot):
    await bot.add_cog(DevCommands(bot))
