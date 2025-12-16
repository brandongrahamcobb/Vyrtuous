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
from vyrtuous.service.discord_message_service import DiscordMessageService
from vyrtuous.service.member_service import MemberService
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.database import Database
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.temporary_room import TemporaryRoom

class DevCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.member_service = MemberService()

    # DONE
    @app_commands.command(name='alias', description='Set an alias for a Vyrtuous action.')
    @is_owner_developer_predicator()
    @app_commands.describe(
        alias_type='One of: cow, uncow, mute, unmute, ban, unban, flag, unflag, tmute, untmute, role, unrole',
        alias_name='Alias/Pseudonym',
        channel='Tag a channel or include its snowflake ID',
        role='Role ID (only for role/unrole)'
    )
    async def create_alias_app_command(
        self,
        interaction: discord.Interaction,
        alias_type: Optional[str] = None,
        alias_name: Optional[str] = None,
        channel: Optional[str] = None,
        role: Optional[str] = None
    ):
        valid_types = {'cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute', 'role', 'unrole'}
        if not alias_type or alias_type.lower() not in valid_types:
            return await interaction.response.send_message(content=f'\U0001F6AB Invalid alias type. Must be one of: `{"`, `".join(valid_types)}`')
        alias_type = alias_type.lower()
        if not alias_name or not alias_name.strip():
            return await interaction.response.send_message(content='\U0001F6AB Alias name cannot be empty.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        old_aliases = await Alias.fetch_command_aliases_by_channel(channel_obj)
        if old_aliases:
            for old_alias in old_aliases:
                if old_alias.alias_name == alias_name:
                    return await interaction.response.send_message(content=f'\U0001F6AB Alias `{alias_name}` ({alias_type}) already exists and is set to  <@{old_alias.channel_id}>.')
        if alias_type in ('role', 'unrole') and not role:
            return await interaction.response.send_message(content='\U0001F6AB Role ID is required for role/unrole aliases.')
        if role:
            role_id = int(role.replace('<@&','').replace('>',''))
        else:
            role_id = None
        alias_obj = Alias(interaction.guild.id, channel_obj.id, alias_type, alias_name, role_id)
        await alias_obj.insert_into_command_aliases()
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_alias', None, interaction.user.id, interaction.guild.id, channel_obj.id, f'Created an alias: {alias_name}')
        if alias_type in ('role','unrole'):
            mention = interaction.guild.get_role(role_id).mention
        else:
            mention = channel_obj.mention
        await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Alias `{alias_name}` ({alias_type}) set to {mention}.', allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @commands.command(name='alias', help='Set an alias for a cow, uncow, mute, unmute, ban, unban, flag, unflag, tmute, untmute, role, or unrole action.')
#    @is_owner_developer_predicator()
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        alias_type: Optional[str] = commands.parameter(default=None, description='One of: `cow`, `uncow`, `mute`, `unmute`, `ban`, `unban`, `flag`, `unflag`, `tmute`, `untmute`, `role`, `unrole`'),
        alias_name: Optional[str] = commands.parameter(default=None, description='Alias/Pseudonym'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        *,
        role: Optional[str] = commands.parameter(default=None, description='Role ID (only for role/unrole)')
    ) -> None:
        valid_types = {'cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'cow', 'uncow', 'tmute', 'untmute', 'role', 'unrole'}
        if not alias_type or alias_type.lower() not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid alias type. Must be one of: `{"`, `".join(valid_types)}`')
        alias_type = alias_type.lower()
        if not alias_name or not alias_name.strip():
            return await self.handler.send_message(ctx, content='\U0001F6AB Alias name cannot be empty.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        old_aliases = await Alias.fetch_command_aliases_by_channel(channel_obj)
        if old_aliases:
            for old_alias in old_aliases:
                if old_alias.alias_name == alias_name:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB Alias `{old_alias.alias_name}` ({old_alias.alias_type}) already exists and is set to <@{old_alias.channel_id}>.')
        if alias_type in ('role', 'unrole') and not role:
            return await self.handler.send_message(ctx, content='\U0001F6AB Role ID is required for role/unrole aliases.')
        if role:
            role_id = int(role.replace('<@&','').replace('>',''))
        else:
            role_id = None
        alias_obj = Alias(ctx.guild.id, channel_obj.id, alias_type, alias_name, role_id)
        await alias_obj.insert_into_command_aliases()
        async with ctx.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_alias', None, ctx.author.id, ctx.guild.id, channel_obj.id, f'Created an alias: {alias_name}')
        if alias_type in ('role','unrole'):
            mention = ctx.guild.get_role(role_id).mention
        else:
            mention = channel_obj.mention
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Alias `{alias_name}` ({alias_type}) set to {mention}.', allowed_mentions=discord.AllowedMentions.none())
    
    # TODO
    @app_commands.command(name='backup', description='Creates a backup of the database and uploads it.')
    @is_owner_developer_predicator()
    async def app_backup(
        self,
        interaction: discord.Interaction
    ):
        try:
            backup = Database(directory='/app/backups')
            backup.create_backup_directory()
            backup_file = backup.execute_backup()
            if backup_file:
                await interaction.response.send_message(file=discord.File(backup.file_name))
            else:
                await interaction.response.send_message(content=f'\U0001F6AB Failed to create backup.')
        except Exception as e:
            logger.warning(f'Database error occurred: {e}')
          
    # TODO
    @commands.command(name='backup', help='Creates a backup of the database and uploads it.')
    @is_owner_developer_predicator()
    async def backup(
        self,
        ctx: commands.Context
    ):
        try:
            backup = Database()
            backup.create_backup_directory()
            backup_file = backup.execute_backup()
            if backup_file:
                await self.handler.send_message(ctx, file=discord.File(backup.file_name))
            else:
                await self.handler.send_message(ctx, content=f'\U0001F6AB Failed to create backup.')
        except Exception as e:
            logger.warning(f'Database error occurred: {e}')
 
    # DONE
    @app_commands.command(name='clear', description='Removes a specific channel ID from all users, including temp-room associations.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_predicator()
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $1::bigint),
                    moderator_channel_ids = array_remove(moderator_channel_ids, $1::bigint),
                    updated_at = NOW()
            ''', channel_obj.id)
            await TemporaryRoom.delete_temporary_room_by_channel(channel_obj)
        Alias.delete_all_command_aliases_by_channel(channel_obj)
        await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Removed channel ID `{channel_obj.id}` from all users\' coordinator and moderator access and deleted all associated records.', allowed_mentions=discord.AllowedMentions.none())

    # DONE
    @commands.command(name='clear', help='Removes a specific channel ID from all users, including temp-room associations and all related records.')
    @is_owner_developer_predicator()
    async def clear_channel_access_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE users
                SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $1::bigint),
                    moderator_channel_ids = array_remove(moderator_channel_ids, $1::bigint),
                    updated_at = NOW()
            ''', channel_obj.id)
            await TemporaryRoom.delete_temporary_room_by_channel(channel_obj)
        Alias.delete_all_command_aliases_by_channel(channel_obj)
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Removed channel {channel_obj.mention} from all users and deleted all associated records.', allowed_mentions=discord.AllowedMentions.none())
    
    @app_commands.command(name='load', description='Loads a cog by name "vyrtuous.cog.<cog_name>.')
    @is_owner_developer_predicator()
    async def load_app_command(self, interaction: discord.Interaction, module: str):
        try:
            if interaction:
                await interaction.response.defer(ephemeral=True)
            await interaction.client.load_extension(module)
        except commands.ExtensionError as e:
            await interaction.response.send_message(f'{e.__class__.__name__}: {e}')
        else:
            await interaction.response.send_message('\N{OK HAND SIGN}')
            
    @commands.command(name='load', help='Loads a cog by name "vyrtuous.cog.<cog_name>."')
    @is_owner_developer_predicator()
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')
    # DONE
    @app_commands.command(name='migrate', description='Migrate a temporary room to a new channel.')
    @app_commands.describe(old_name='Old temporary room name', new_channel='New channel to migrate to')
    @is_owner_developer_predicator()
    async def migrate_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        old_name: str,
        new_channel: discord.abc.GuildChannel
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        async with self.bot.db_pool.acquire() as conn:
            room = await TemporaryRoom.fetch_temporary_room_by_guild_and_room_name(interaction.guild, old_name)
            if not room:
                return await interaction.response.send_message(content=f'\U0001F6AB No temporary room named `{old_name}` found.')
            channel_obj = await self.channel_service.resolve_channel(interaction, new_channel.id)
            if channel_obj.type != discord.ChannelType.voice:
                return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
            is_owner = room.room_owner == interaction.user
            await room.update_temporary_room_name_and_room_snowflake(channel_obj, new_channel.name)
            aliases = await Alias.fetch_command_aliases_by_channel_id(interaction.guild.id, room.room_snowflake)
            if aliases:
                for alias_obj in aliases:
                    await alias_obj.update_command_aliases_with_channel(channel_obj)
            tables = [
                'active_bans','active_text_mutes','active_voice_mutes',
                'active_stages','stage_coordinators','active_caps'
            ]
            for table in tables:
                await conn.execute(
                    f'UPDATE {table} SET room_name=$3, channel_id=$4 WHERE guild_id=$1 AND room_name=$2',
                    interaction.guild.id, old_name, channel_obj.name, channel_obj.id
                )
            await conn.execute(
                'UPDATE users SET coordinator_channel_ids=array_replace(coordinator_channel_ids, $1::bigint, $2::bigint) WHERE $1=ANY(coordinator_channel_ids)',
                room.room_snowflake, channel_obj.id
            )
            await conn.execute(
                'UPDATE users SET moderator_channel_ids=array_replace(moderator_channel_ids, $1::bigint, $2::bigint) WHERE $1=ANY(moderator_channel_ids)',
                room.room_snowflake, channel_obj.id
            )
            return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Temporary room `{old_name}` migrated to {channel_obj.mention} and renamed to `{channel_obj.name}`.')
    
    # DONE
    @commands.command(name='migrate', help='Migrate a temporary room to a new channel by snowflake.')
    @is_owner_developer_predicator()
    async def migrate_temp_room_text_command(
        self,
        ctx,
        old_name: Optional[str] = commands.parameter(default=None, description='Provide a channel name'),
        new_channel: Optional[int] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        async with self.bot.db_pool.acquire() as conn:
            room = await TemporaryRoom.fetch_temporary_room_by_guild_and_room_name(ctx.guild, old_name)
            if not room:
                return await self.handler.send_message(ctx, content=f'No temporary room named `{old_name}` found.')
            channel_obj = await self.channel_service.resolve_channel(ctx, new_channel)
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            is_owner = room.room_owner == ctx.author
            await room.update_temporary_room_name_and_room_snowflake(channel_obj, channel_obj.name)
            aliases = await Alias.fetch_command_aliases_by_channel_id(ctx.guild.id, room.room_snowflake)
            if aliases:
                for alias_obj in aliases:
                    await alias_obj.update_command_aliases_with_channel(channel_obj)
            tables = [
                'active_bans','active_text_mutes','active_voice_mutes',
                'active_stages','stage_coordinators','active_caps'
            ]
            for table in tables:
                await conn.execute(
                    f'UPDATE {table} SET room_name=$3, channel_id=$4 WHERE guild_id=$1 AND room_name=$2',
                    ctx.guild.id, old_name, channel_obj.name, channel_obj.id
                )
            await conn.execute(
                'UPDATE users SET coordinator_channel_ids=array_replace(coordinator_channel_ids, $1::bigint, $2::bigint) WHERE $1=ANY(coordinator_channel_ids)',
                room.room_snowflake, channel_obj.id
            )
            await conn.execute(
                'UPDATE users SET moderator_channel_ids=array_replace(moderator_channel_ids, $1::bigint, $2::bigint) WHERE $1=ANY(moderator_channel_ids)',
                room.room_snowflake, channel_obj.id
            )
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Temporary room `{old_name}` migrated to {channel_obj.mention} and renamed to `{channel_obj.name}`.')

    @app_commands.command(name='reload', description='Reloads a cog by name "vyrtuous.cog.<cog_name>".')
    @app_commands.check(at_home)
    @is_owner_developer_predicator()
    async def reload_app_command(self, interaction: discord.Interaction, module: str):
        try:
            if interaction:
                await interaction.response.defer(ephemeral=True)
            await interaction.client.reload_extension(module)
        except commands.ExtensionError as e:
            await interaction.response.send_message(f'{e.__class__.__name__}: {e}')
        else:
            await interaction.response.send_message('\N{OK HAND SIGN}')
            
    @commands.command(name='reload', help='Reloads a cog by name "vyrtuous.cog.<cog_name>".')
    @is_owner_developer_predicator()
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @app_commands.command(name='sync', description="Syncs commands to the tree '~', globally '*', clear '^' or general sync.")
    @is_owner_developer_predicator()
    async def sync_app_command(self, interaction: discord.Interaction, spec: Optional[Literal['~', '*', '^']] = None) -> None:
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
            await interaction.followup.send(f"Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}")
            return
        ret = 0
        for guild in guilds:
            try:
                await interaction.client.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        await interaction.followup.send(f'Synced the tree to {ret}/{len(guilds)}.')
        
    @commands.command(name='sync', help="Syncs commands to the tree '~', globally '*', clear '^' or general sync.")
    @is_owner_developer_predicator()
    async def sync_text_command(self, ctx: commands.Context, guilds: commands.Greedy[discord.Object], spec: Optional[Literal['~', '*', '^']] = None) -> None:
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
            await ctx.send(f"Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}")
            return
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        await ctx.send(f'Synced the tree to {ret}/{len(guilds)}.')

    # DONE
    @app_commands.command(name='temp', description='Toggle a temporary room and assign an owner.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID', owner='Tag a member or include their snowflake ID')
    @is_owner_developer_predicator()
    async def toggle_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None,
        owner: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if not isinstance(channel_obj, discord.VoiceChannel):
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        room = await TemporaryRoom.fetch_temporary_room_by_channel(channel_obj)
        action = None
        if room:
            async with self.bot.db_pool.acquire() as conn:
                if room.room_owner:
                    await conn.execute(
                        'UPDATE users SET moderator_channel_ids=array_remove(moderator_channel_ids,$1), updated_at=NOW() '
                        'WHERE discord_snowflake=$2',
                        channel_obj.id,
                        room.room_owner.id
                    )
                await TemporaryRoom.delete_temporary_room_by_channel(channel_obj)
                await Alias.delete_all_command_aliases_by_channel(channel_obj)
            action = 'removed'
        else:
            member_obj = await self.member_service.resolve_member(interaction, owner)
            temporary_room = TemporaryRoom(interaction.guild, channel_obj.id, member_obj)
            async with self.bot.db_pool.acquire() as conn:
                await temporary_room.insert_into_temporary_rooms()
                await conn.execute(
                    'UPDATE users SET moderator_channel_ids = CASE WHEN NOT $1=ANY(moderator_channel_ids) '
                    'THEN array_append(moderator_channel_ids,$1) ELSE moderator_channel_ids END, updated_at=NOW() '
                    'WHERE discord_snowflake=$2',
                    channel_obj.id,
                    member_obj.id
                )
            action = f'created and owned by {member_obj.mention}'
        await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Temporary room {channel_obj.mention} has been {action}.', allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @commands.command(name='temp', help='Toggle a temporary room and assign an owner.')
    @is_owner_developer_predicator()
    async def toggle_temp_room_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        owner: Optional[str] = commands.parameter(default=None, description='Tag a member or include their Discord ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if not isinstance(channel_obj, discord.VoiceChannel):
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        room = await TemporaryRoom.fetch_temporary_room_by_channel(channel_obj)
        action = None
        if room:
            async with ctx.bot.db_pool.acquire() as conn:
                if room.room_owner:
                    await conn.execute(
                        'UPDATE users SET moderator_channel_ids=array_remove(moderator_channel_ids,$1), updated_at=NOW() '
                        'WHERE discord_snowflake=$2',
                        channel_obj.id,
                        room.room_owner.id
                    )
                await TemporaryRoom.delete_temporary_room_by_channel(channel_obj)
                await Alias.delete_all_command_aliases_by_channel(channel_obj)
            action = 'removed'
        else:
            member_obj = await self.member_service.resolve_member(ctx, owner)
            temporary_room = TemporaryRoom(ctx.guild, channel_obj.id, member_obj)
            async with ctx.bot.db_pool.acquire() as conn:
                await temporary_room.insert_into_temporary_rooms()
                await conn.execute(
                    'UPDATE users SET moderator_channel_ids = CASE WHEN NOT $1=ANY(moderator_channel_ids) '
                    'THEN array_append(moderator_channel_ids,$1) ELSE moderator_channel_ids END, updated_at=NOW() '
                    'WHERE discord_snowflake=$2',
                    channel_obj.id,
                    member_obj.id
                )
            action = f'created and owned by {member_obj.mention}'
        await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Temporary room {channel_obj.mention} has been {action}.', allowed_mentions=discord.AllowedMentions.none())

    # DONE
    @app_commands.command(name='trole', description='Marks a role as administrator and syncs all members.')
    @is_owner_developer_predicator()
    async def grant_administrator_to_role_app_command(self, interaction: discord.Interaction, role: Optional[str]):
        guild = interaction.guild
        if not guild:
            return await interaction.response.send_message(content='\U0001F6AB This command must be used in a server.')
        role_id = int(role.replace('<@&','').replace('>',''))
        role_obj = guild.get_role(role_id)
        async with self.bot.db_pool.acquire() as conn:
            for member in role_obj.members:
                await conn.execute(
                    '''
                    INSERT INTO users (discord_snowflake, administrator_role_ids, administrator_guild_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[], ARRAY[$3]::BIGINT[])
                    ON CONFLICT (discord_snowflake)
                    DO UPDATE SET
                        administrator_role_ids = ARRAY(SELECT DISTINCT unnest(users.administrator_role_ids || $2)),
                        administrator_guild_ids = ARRAY(SELECT DISTINCT unnest(users.administrator_guild_ids || $3))
                    ''',
                    member.id,
                    role_obj.id,
                    guild.id
                )
        await interaction.response.send_message(f'Administrator role `{role_obj.name}` synced.')
    
    @commands.command(name='trole', help='Marks a role as administrator and syncs all members.')
    @is_owner_developer_predicator()
    async def grant_administrator_to_role_text_command(self, ctx: commands.Context, role: Optional[str]):
        guild = ctx.guild
        role_id = int(role.replace('<@&','').replace('>',''))
        role_obj = guild.get_role(role_id)
        async with self.bot.db_pool.acquire() as conn:
            for member in role_obj.members:
                await conn.execute(
                    '''
                    INSERT INTO users (discord_snowflake, administrator_role_ids, administrator_guild_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[], ARRAY[$3]::BIGINT[])
                    ON CONFLICT (discord_snowflake)
                    DO UPDATE SET
                        administrator_role_ids = ARRAY(SELECT DISTINCT unnest(users.administrator_role_ids || $2)),
                        administrator_guild_ids = ARRAY(SELECT DISTINCT unnest(users.administrator_guild_ids || $3))
                    ''',
                    member.id,
                    role_obj.id,
                    guild.id
                )
        await self.handler.send_message(ctx, content=f'Administrator role `{role_obj.name}` synced.')

    # DONE
    @app_commands.command(name='unload', description='Unloads a cog by name "vyrtuous.cog.<cog_name>".')
    @is_owner_developer_predicator()
    async def unload_app_command(self, interaction: discord.Interaction, module: str):
        try:
            await interaction.client.unload_extension(module)
        except commands.ExtensionError as e:
            await interaction.response.send_message(f'{e.__class__.__name__}: {e}')
        else:
            await interaction.response.send_message('\N{OK HAND SIGN}')

    # DONE
    @commands.command(name='unload', help='Unload a cog by name "vyrtuous.cog.<cog_name>".')
    @is_owner_developer_predicator()
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.unload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')
                
    # DONE
    @app_commands.command(name='xalias', description='Deletes an alias.')
    @is_owner_developer_predicator()
    @app_commands.describe(alias_name='Include an alias name')
    async def delete_alias_app_command(
        self,
        interaction: discord.Interaction,
        alias_name: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        if not alias_name or not alias_name.strip():
            return await interaction.response.send_message(content='\U0001F6AB `alias_name` cannot be empty.')
        aliases = await Alias.fetch_command_aliases_by_guild(interaction.guild)
        if not aliases:
            return await interaction.response.send_message(content=f'\U0001F6AB No aliases found.')
        found = False
        for alias in aliases:
            if alias.alias_name:
                await Alias.delete_command_alias_by_guild_and_alias_name(interaction.guild, alias.alias_name)
                found = True
                break
        if not found:
            return await interaction.response.send_message(content=f'\U0001F6AB Alias `{alias_name}` not found.')
        channel_obj = await self.channel_service.resolve_channel(interaction, alias.channel_id)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1,$2,$3,$4,$5,$6)
                ''', 'delete_alias', None, interaction.user.id, interaction.guild.id, channel_obj.id, f'Deleted alias {alias_name}')
        return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Deleted alias `{alias_name}` from `{alias.alias_type}`.')

    # DONE
    @commands.command(name='xalias', help='Deletes an alias.')
    @is_owner_developer_predicator()
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: Optional[str] = commands.parameter(default=None, description='Include an alias name')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        if not alias_name or not alias_name.strip():
            return await self.handler.send_message(ctx, content='\U0001F6AB `alias_name` cannot be empty.')
        aliases = await Alias.fetch_command_aliases_by_guild(ctx.guild)
        if not aliases:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB No aliases found.')
        found = False
        for alias in aliases:
            if alias.alias_name:
                await Alias.delete_command_alias_by_guild_and_alias_name(ctx.guild, alias.alias_name)
                found = True
                break
        if not found:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Alias `{alias_name}` not found.')
        channel_obj = await self.channel_service.resolve_channel(ctx, alias.channel_id)
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs(action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES($1,$2,$3,$4,$5,$6)
            ''', 'delete_alias', None, ctx.author.id, ctx.guild.id, channel_obj.id, f'Deleted alias {alias_name}')
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Deleted alias `{alias_name}` from `{alias.alias_type}`.')

    @app_commands.command(name='xtrole', description='Revokes a role from administrator and updates all members.')
    @is_owner_developer_predicator()
    async def revoke_administrator_to_role_app_command(self, interaction: discord.Interaction, role: Optional[str]):
        guild = interaction.guild
        if not guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        role_id = int(role.replace('<@&','').replace('>',''))
        role_obj = guild.get_role(role_id)
        async with self.bot.db_pool.acquire() as conn:
            for member in role_obj.members:
                row = await conn.fetchrow(
                    'SELECT administrator_role_ids, administrator_guild_ids FROM users WHERE discord_snowflake=$1',
                    member.id
                )
                if not row:
                    continue
                role_ids = set(row['administrator_role_ids'] or [])
                guild_ids = set(row['administrator_guild_ids'] or [])
                role_ids.discard(role_obj.id)
                if not (role_ids & {r.id for r in member.roles}):
                    guild_ids.discard(guild.id)
                await conn.execute(
                    'UPDATE users SET administrator_role_ids=$2, administrator_guild_ids=$3 WHERE discord_snowflake=$1',
                    member.id,
                    list(role_ids),
                    list(guild_ids)
                )
        await interaction.response.send_message(f'Administrator role `{role_obj.name}` removed.')
     
    @commands.command(name='xtrole', help='Revokes a role from administrator and updates all members.')
    @is_owner_developer_predicator()
    async def revoke_administrator_to_role_text_command(self, ctx: commands.Context, role: Optional[str]):
        guild = ctx.guild
        role_id = int(role.replace('<@&','').replace('>',''))
        role_obj = guild.get_role(role_id)
        async with self.bot.db_pool.acquire() as conn:
            for member in role_obj.members:
                row = await conn.fetchrow(
                    'SELECT administrator_role_ids, administrator_guild_ids FROM users WHERE discord_snowflake=$1',
                    member.id
                )
                if not row:
                    continue
                role_ids = set(row['administrator_role_ids'] or [])
                guild_ids = set(row['administrator_guild_ids'] or [])
                role_ids.discard(role_obj.id)
                if not (role_ids & {r.id for r in member.roles}):
                    guild_ids.discard(guild.id)
                await conn.execute(
                    'UPDATE users SET administrator_role_ids=$2, administrator_guild_ids=$3 WHERE discord_snowflake=$1',
                    member.id,
                    list(role_ids),
                    list(guild_ids)
                )
        await self.handler.send_message(ctx, content=f'Administrator role `{role_obj.name}` removed.')         
async def setup(bot: DiscordBot):
    await bot.add_cog(DevCommands(bot))
