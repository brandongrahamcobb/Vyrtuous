''' admin_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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
from vyrtuous.inc.helpers import *
from vyrtuous.service.check_service import *
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.check_service import *
from vyrtuous.service.member_service import MemberService
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.duration import Duration
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.statistics import Statistics
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.service.discord_message_service import AppPaginator, DiscordMessageService, Paginator

class AdminCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.member_service = MemberService()
        
    @app_commands.command(name='cap', description='Set a duration limit for bans, mutes and text mutes.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(
        channel='Tag a channel or include its snowflake ID',
        moderation_type='One of: `mute`, `ban`, `tmute`',
        duration='(+|-)duration(m|h|d), 0=permanent, default=24h'
    )
    async def cap_app_command(
        self,
        interaction: discord.Interaction,
        channel: str = None,
        moderation_type: str = None,
        duration: str = '24h'
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.') 
        valid_types = {'mute', 'ban', 'tmute'}
        if moderation_type not in valid_types:
            return await interaction.response.send_message(content=f'\U0001F6AB Invalid moderation type. Must be one of: {", ".join(valid_types)}')
        duration_obj = Duration()
        duration_obj.load_from_combined_duration_str(duration)
        expires_at = duration_obj.output_datetime()
        duration_display = duration_obj.output_display()
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4',
                interaction.guild.id, channel_obj.id, moderation_type, channel_obj.name
            )
            original_duration = row['duration'] if row else None
            await conn.execute('''
                INSERT INTO active_caps (guild_id, channel_id, moderation_type, duration, room_name)
                VALUES ($1, $2, $3, $4, $5)
                ON CONFLICT (guild_id, channel_id, moderation_type, room_name)
                DO UPDATE SET duration = EXCLUDED.duration
            ''', interaction.guild.id, channel_obj.id, moderation_type, duration, channel_obj.name)
        if original_duration:
            msg = f'{self.emoji.get_random_emoji()} Cap changed on {channel_obj.mention} for {moderation_type} from {original_duration} to {duration_display}.'
        else:
            msg = f'{self.emoji.get_random_emoji()} Cap set on {channel_obj.mention} for {moderation_type} for {duration_display}.'
        await interaction.response.send_message(content=msg, allowed_mentions=discord.AllowedMentions.none())
    
    # DONE
    @commands.command(name='cap', help='Set a duration limit for bans, mutes and text mutes.')
    @is_owner_developer_administrator_predicator()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        moderation_type: Optional[str] = commands.parameter(default=None, description='One of: `mute`, `ban`, `tmute`'),
        *,
        duration: Optional[str] = commands.parameter(default='24h', description='Options: (+|-)duration(m|h|d) \n 0 - permanent / 24h - default')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        valid_types = {'mute', 'ban', 'tmute'}
        if moderation_type not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid moderation type. Must be one of: {", ".join(valid_types)}')
        duration_obj = Duration()
        duration_obj.load_from_combined_duration_str(duration)
        expires_at = duration_obj.output_datetime()
        duration_display = duration_obj.output_display()
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4',
                ctx.guild.id, channel_obj.id, moderation_type, channel_obj.name
            )
            original_duration = row['duration'] if row else None
            await conn.execute('''
                INSERT INTO active_caps (guild_id, channel_id, moderation_type, duration, room_name)
                VALUES ($1, $2, $3, $4, $5)
                ON CONFLICT (guild_id, channel_id, moderation_type, room_name)
                DO UPDATE SET duration = EXCLUDED.duration
            ''', ctx.guild.id, channel_obj.id, moderation_type, duration, channel_obj.name)
        if original_duration:
            msg = f'{self.emoji.get_random_emoji()} Cap changed on {channel_obj.mention} for {moderation_type} from {original_duration} to {duration_display}.'
        else:
            msg = f'{self.emoji.get_random_emoji()} Cap set on {channel_obj.mention} for {moderation_type} for {duration_display}.'
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
    
    # DONE
    @app_commands.command(name='chown', description='Change the owner of a temporary room.')
    @app_commands.describe(member='Tag a user or provide their snowflake ID', channel='Tag a channel or provide it\'s snowflake ID')
    @is_owner_developer_administrator_predicator()
    async def change_temp_room_owner_app_command(
        self,
        interaction,
        channel: Optional[str],
        member: Optional[str]
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(interaction, member)
        if not member_obj:
            return await interaction.response.send_message(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        room = await TemporaryRoom.fetch_temporary_room_by_channel(channel_obj)
        if not room:
            return await interaction.response.send_message(content=f'{channel_obj.mention} is not a temporary room.')
        await room.update_temporary_room_owner_snowflake(member_obj)
        return await interaction.response.send_message(content=f'Ownership of {channel_obj.mention} has been transferred to {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            
    # DONE
    @commands.command(name='chown', help='Change the owner of a temporary room.')
    @is_owner_developer_administrator_predicator()
    async def change_temp_room_owner_text_command(
        self,
        ctx,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or provide it\'s snowflake ID'),
        member: Optional[str] = commands.parameter(default=None, description='Tag a user or provide their snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        room = await TemporaryRoom.fetch_temporary_room_by_channel(channel_obj)
        if not room:
            return await self.handler.send_message(ctx, content=f'{channel_obj.mention} is not a temporary room.')
        await room.update_temporary_room_owner_snowflake(member_obj)
        return await self.handler.send_message(ctx, content=f'Ownership of {channel_obj.mention} has been transferred to {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @app_commands.command(name='coord', description='Grants/revokes coordinator access for a specific voice channel.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(member='Tag a member or include their snowflake ID', channel='Tag a channel or include its snowflake ID')
    async def create_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(interaction, member)
        if not member_obj:
            return await self.send(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        if member_obj.bot:
            return await interaction.response.send_message(content='\U0001F6AB You cannot make the bot a coordinator.')
        success = await has_equal_or_higher_role(interaction, member_obj, channel_obj)
        if not success:
            return await interaction.response.send_message(content=f"\U0001F6AB You are not allowed to toggle {member_obj.mention}'s role as a coordinator because they are a higher/or equivalent role than you in {channel_obj.mention}.", allowed_mentions=discord.AllowedMentions.none())
        async with self.bot.db_pool.acquire() as conn:
            action = None
            row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            current_channel_ids = row.get('coordinator_channel_ids', []) if row else []
            if channel_obj.id in current_channel_ids:
                await conn.execute('UPDATE users SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $2), updated_at = NOW() WHERE discord_snowflake = $1', member_obj.id, channel_obj.id)
                action = 'revoked'
            else:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, coordinator_channel_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET coordinator_channel_ids = (
                        SELECT ARRAY(SELECT DISTINCT unnest(
                            COALESCE(users.coordinator_channel_ids, ARRAY[]::BIGINT[]) || ARRAY[$2]::BIGINT[]
                        ))
                    ),
                    updated_at = NOW()
                ''', member_obj.id, channel_obj.id)
                action = 'granted'
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_coordinator', member_obj.id, interaction.user.id, interaction.guild.id, channel_obj.id, f'Coordinator access {action}')
        return await interaction.response.send_message(content=f"{self.emoji.get_random_emoji()} {member_obj.mention}'s coordinator access has been {action} in {channel_obj.mention}.", allowed_mentions=discord.AllowedMentions.none())

    # DONE
    @commands.command(name='coord', help='Grants/revokes coordinator access for a specific voice channel.')
    @is_owner_developer_administrator_predicator()
    async def create_coordinator_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        if member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a coordinator.')
        success = await has_equal_or_higher_role(ctx.message, member_obj, channel_obj)
        if not success:
            return await self.handler.send_message(ctx, content=f"\U0001F6AB You are not toggle {member_obj.mention}'s coordinator role because they are a higher/or equivalent role than you in {channel_obj.mention}.", allowed_mentions=discord.AllowedMentions.none())
        async with self.bot.db_pool.acquire() as conn:
            action = None
            row = await conn.fetchrow('SELECT coordinator_channel_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            current_channel_ids = row.get('coordinator_channel_ids', []) if row else []
            if not current_channel_ids:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, coordinator_channel_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET coordinator_channel_ids = (
                        SELECT ARRAY(SELECT DISTINCT unnest(
                            COALESCE(users.coordinator_channel_ids, ARRAY[]::BIGINT[]) || ARRAY[$2]::BIGINT[]
                        ))
                    ),
                    updated_at = NOW()
                ''', member_obj.id, channel_obj.id)
                action = 'granted'
            elif channel_obj.id in current_channel_ids:
                await conn.execute('UPDATE users SET coordinator_channel_ids = array_remove(coordinator_channel_ids, $2), updated_at = NOW() WHERE discord_snowflake = $1', member_obj.id, channel_obj.id)
                action = 'revoked'
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_coordinator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, f'Coordinator access {action}')
        return await self.handler.send_message(ctx, content=f"{self.emoji.get_random_emoji()} {member_obj.mention}'s coordinator access has been {action} in {channel_obj.mention}.", allowed_mentions=discord.AllowedMentions.none())
        
   # DONE
    @app_commands.command(name='cstage', description='Create a stage in the current or specified channel.')
    @app_commands.describe(channel='Tag a voice/stage channel', duration='Duration of the stage (e.g., 1h, 30m)')
    @is_owner_developer_administrator_predicator()
    async def stage_create_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None,
        duration: str = '1'
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content="This command must be used in a server.")
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        duration_obj = Duration()
        duration_obj.load_from_combined_duration_str(duration)
        expires_at = duration_obj.output_datetime()
        duration_display = duration_obj.output_display()
        skipped, muted, failed = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (guild_id, channel_id, initiator_id, expires_at, room_name)
                VALUES ($1, $2, $3, $4, $5)
                ON CONFLICT (guild_id, channel_id, room_name)
                DO UPDATE SET expires_at=EXCLUDED.expires_at, initiator_id=EXCLUDED.initiator_id
            ''', interaction.guild.id, channel_obj.id, interaction.user.id, expires_at, channel_obj.name)
            await conn.execute('''
                INSERT INTO stage_coordinators (guild_id, channel_id, discord_snowflake, room_name)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT DO NOTHING
            ''', interaction.guild.id, channel_obj.id, interaction.user.id, channel_obj.name)
            for user in channel_obj.members:
                if await member_is_owner(user) or await member_is_developer(user) or await member_is_administrator(user) \
                   or await member_is_coordinator(user, channel_obj) or await member_is_moderator(user, channel_obj) \
                   or user.id == interaction.user.id:
                    skipped.append(user)
                    continue
                try:
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason, target, room_name)
                        VALUES ($1,$2,$3,$4,$5,'room',$6)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                        DO UPDATE SET expires_at=EXCLUDED.expires_at
                    ''', interaction.guild.id, user.id, channel_obj.id, expires_at, 'Stage active', channel_obj.name)
                    if user.voice and user.voice.channel.id == channel_obj.id:
                        await user.edit(mute=True)
                    muted.append(user)
                except Exception as e:
                    logger.warning(f'Failed to mute {user.name}: {e}')
                    failed.append(user)
        msg = f'{self.emoji.get_random_emoji()} Stage created for {duration_display} in {channel_obj.mention}.\nMuted {len(muted)} user(s). Skipped {len(skipped)}.'
        if failed:
            msg += f'\U0001F6AB Failed: {len(failed)}.'
        await interaction.response.send_message(content=msg, allowed_mentions=discord.AllowedMentions.none())
    
    # DONE
    @commands.command(name='cstage', help='Create a stage in the current or specified channel.')
    @is_owner_developer_administrator_predicator()
    async def stage_create_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        *,
        duration: str = '1'
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        duration_obj = Duration()
        duration_obj.load_from_combined_duration_str(duration)
        expires_at = duration_obj.output_datetime()
        duration_display = duration_obj.output_display()
        skipped, muted, failed = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_stages (guild_id, channel_id, initiator_id, expires_at, room_name)
                VALUES ($1, $2, $3, $4, $5)
                ON CONFLICT (guild_id, channel_id, room_name)
                DO UPDATE SET expires_at=EXCLUDED.expires_at, initiator_id=EXCLUDED.initiator_id
            ''', ctx.guild.id, channel_obj.id, ctx.author.id, expires_at, channel_obj.name)
            await conn.execute('''
                INSERT INTO stage_coordinators (guild_id, channel_id, discord_snowflake, room_name)
                VALUES ($1, $2, $3, $4)
                ON CONFLICT DO NOTHING
            ''', ctx.guild.id, channel_obj.id, ctx.author.id, channel_obj.name)
            for user in channel_obj.members:
                if await member_is_owner(user) or await member_is_developer(user) or await member_is_administrator(user) \
                   or await member_is_coordinator(user, channel_obj) or await member_is_moderator(user, channel_obj) \
                   or user.id == ctx.author.id:
                    skipped.append(user)
                    continue
                try:
                    await conn.execute('''
                        INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason, target, room_name)
                        VALUES ($1,$2,$3,$4,$5,'room',$6)
                        ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                        DO UPDATE SET expires_at=EXCLUDED.expires_at
                    ''', ctx.guild.id, user.id, channel_obj.id, expires_at, 'Stage active', channel_obj.name)
                    if user.voice and user.voice.channel.id == channel_obj.id:
                        await user.edit(mute=True)
                    muted.append(user)
                except Exception as e:
                    logger.warning(f'Failed to mute {user.name}: {e}')
                    failed.append(user)
        msg = f'{self.emoji.get_random_emoji()} Stage created for {duration_display} in {channel_obj.mention}.\nMuted {len(muted)} user(s). Skipped {len(skipped)}.'
        if failed:
            msg += f'\n\U000026A0\U0000FE0F Failed: {len(failed)}.'
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @app_commands.command(name='log', description='Toggle logging for a channel on or off.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    async def toggle_log_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid text channel.')
        current_channels = Statistics.get_statistic_channels().setdefault(interaction.guild.id, [])
        async with self.bot.db_pool.acquire() as conn:
            existing = await conn.fetchrow(
                'SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2;',
                interaction.guild.id, channel_obj.id
            )
            if existing:
                new_status = not existing['enabled']
                await conn.execute(
                    'UPDATE statistic_channels SET enabled=$1 WHERE guild_id=$2 AND channel_id=$3;',
                    new_status, interaction.guild.id, channel_obj.id
                )
                msg = f'{self.emoji.get_random_emoji()} Logging for {channel_obj.mention} toggled {"on" if new_status else "off"}.'
            else:
                msg = f'\U0001F6AB No logging for {channel_obj.mention} exists'
        await Statistics.load_channels()
        await interaction.response.send_message(content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @commands.command(name='log', help='Toggle logging for a channel on or off.')
    @is_owner_developer_administrator_predicator()
    async def toggle_log_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid text channel.')
        current_channels = Statistics.get_statistic_channels().setdefault(ctx.guild.id, [])
        async with self.bot.db_pool.acquire() as conn:
            existing = await conn.fetchrow(
                'SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2;',
                ctx.guild.id, channel_obj.id
            )
            if existing:
                new_status = not existing['enabled']
                await conn.execute(
                    'UPDATE statistic_channels SET enabled=$1 WHERE guild_id=$2 AND channel_id=$3;',
                    new_status, ctx.guild.id, channel_obj.id
                )
                msg = f'{self.emoji.get_random_emoji()} Logging for {channel_obj.mention} toggled {"on" if new_status else "off"}.'
            else:
                msg = f'\U0001F6AB No logging for {channel_obj.mention} exists'
        await Statistics.load_channels()
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @app_commands.command(name='logs', description='Lists log channels for this guild or a specified guild ID.')
    @is_owner_developer_administrator_predicator()
    async def list_logs_app_command(
        self,
        interaction: discord.Interaction
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        entries = Statistics.get_statistic_channels().setdefault(interaction.guild.id, [])
        if not entries:
            return await interaction.response.send_message(content=f'\U0001F6AB No log channels configured in {interaction.guild.name}.')
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Log Channels',
            color=discord.Color.blue()
        )
        embed.set_footer(text=f'Guild ID: {interaction.guild.id}')
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1;', interaction.guild.id)
            for row in rows:
                channel_obj = self.bot.get_channel(row['channel_id'])
                mention = channel_obj.mention if channel_obj else f'`{row["channel_id"]}`'
                enabled = row.get('enabled',False)
                log_type = row.get('type') or 'general'
                snowflakes = row.get('snowflakes') or []
                if log_type == 'general':
                    detail = f"Logs all events in {interaction.guild.name}"
                elif log_type == 'channel':
                    detail=f"Logs only events from channels: {', '.join(f'<#{s}>' for s in snowflakes)}"
                elif log_type == 'member':
                    detail=f"Logs only events for members: {', '.join(f'<@{s}>' for s in snowflakes)}"
                else:
                    detail="Unknown filter"
                embed.add_field(name=f"{mention} {'✅' if enabled else '\U0001F6AB'}", value=f"Type: **{log_type}**\n{detail}", inline=False)
        await interaction.response.send_message(embed=embed) 
               
    # DONE
    @commands.command(name='logs', help='Lists log channels for this guild or a specified guild ID.')
    @is_owner_developer_administrator_predicator()
    async def list_logs_text_command(
        self,
        ctx: commands.Context
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        entries = Statistics.get_statistic_channels().setdefault(ctx.guild.id, [])
        if not entries:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB No log channels configured in {ctx.guild.name}.')
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Log Channels',
            color=discord.Color.blue()
        )
        embed.set_footer(text=f'Guild ID: {ctx.guild.id}')
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT * FROM statistic_channels WHERE guild_id=$1;', ctx.guild.id)
            for row in rows:
                channel_obj = self.bot.get_channel(row['channel_id'])
                mention = channel_obj.mention if channel_obj else f'`{row["channel_id"]}`'
                enabled = row.get('enabled', False)
                log_type = row.get('type') or 'general'
                snowflakes = row.get('snowflakes') or []
                if log_type == 'general':
                    detail = f"Logs all events in {ctx.guild.name}"
                elif log_type == 'channel':
                    detail = f"Logs only events from channels: {', '.join(f'<#{s}>' for s in snowflakes)}"
                elif log_type == 'member':
                    detail = f"Logs only events for members: {', '.join(f'<@{s}>' for s in snowflakes)}"
                else:
                    detail = "Unknown filter"
                embed.add_field(
                    name=f"{mention} {'✅' if enabled else '\U0001F6AB'}",
                    value=f"Type: **{log_type}**\n{detail}",
                    inline=False
                )
        await self.handler.send_message(ctx, embed=embed)

    # DONE
    @app_commands.command(name='mlog', description='Create, modify, or delete a log channel.')
    @app_commands.describe(
        channel='Tag a channel or include its snowflake ID',
        action='create | modify | delete',
        log_type='Type of logs: member, channel, general',
        snowflakes='Optional list of member IDs to include in logs'
    )
    @is_owner_developer_administrator_predicator()
    async def modify_log_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None,
        action: Optional[str] = None,
        log_type: Optional[str] = None,
        snowflakes: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        if action not in ['create', 'modify', 'delete']:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid action `create`, `modify` or `delete`.')
        sf = [int(s) for s in snowflakes.split()] if snowflakes else []
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid text channel.')
        current_entries = Statistics.get_statistic_channels().setdefault(interaction.guild.id, [])
        async with self.bot.db_pool.acquire() as conn:
            if action.lower() == 'delete':
                await conn.execute('DELETE FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2;', interaction.guild.id, channel_obj.id)
                current_entries[:] = [e for e in current_entries if e['channel_id'] != channel_obj.id]
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Log channel {channel_obj.mention} deleted.')
            existing = await conn.fetchrow('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2;', interaction.guild.id, channel_obj.id)
            if existing:
                await conn.execute('UPDATE statistic_channels SET type=$1, snowflakes=$2, enabled=TRUE WHERE guild_id=$3 AND channel_id=$4;', log_type, sf if sf else None, interaction.guild.id, channel_obj.id)
                for e in current_entries:
                    if e['channel_id'] == channel_obj.id:
                        e.update({'type': log_type, 'snowflakes': sf if sf else None, 'enabled': True})
                msg = f'{self.emoji.get_random_emoji()} Log channel {channel_obj.mention} updated with type `{log_type or "general"}`.'
            else:
                await conn.execute('INSERT INTO statistic_channels (guild_id, channel_id, type, snowflakes, enabled) VALUES ($1, $2, $3, $4, TRUE);', interaction.guild.id, channel_obj.id, log_type, sf if sf else None)
                current_entries.append({'guild_id': interaction.guild.id, 'channel_id': channel_obj.id, 'type': log_type, 'snowflakes': sf if sf else None, 'enabled': True})
                msg = f'{self.emoji.get_random_emoji()} Log channel {channel_obj.mention} created with type `{log_type or "general"}`.'
        await Statistics.load_channels()
        await interaction.response.send_message(content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @commands.command(name='mlog', help='Create, modify, or delete a log channel.')
    @is_owner_developer_administrator_predicator()
    async def modify_log_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        action: Optional[str] = commands.parameter(default=None, description='create | modify | delete'),
        log_type: Optional[str] = commands.parameter(default=None, description='Type of logs: member, channel, general'),
        *snowflakes: Optional[int]
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        if action not in ['create', 'modify', 'delete']:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid action `create`, `modify` or `delete`.')
        sf = [int(s) for s in snowflakes] if snowflakes else []
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.text:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid text channel.')
        current_entries = Statistics.get_statistic_channels().setdefault(ctx.guild.id, [])
        async with self.bot.db_pool.acquire() as conn:
            if action.lower() == 'delete':
                await conn.execute('DELETE FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2;', ctx.guild.id, channel_obj.id)
                current_entries[:] = [e for e in current_entries if e['channel_id'] != channel_obj.id]
                await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Log channel {channel_obj.mention} deleted.')
                return
            existing = await conn.fetchrow('SELECT * FROM statistic_channels WHERE guild_id=$1 AND channel_id=$2;', ctx.guild.id, channel_obj.id)
            if existing:
                await conn.execute('UPDATE statistic_channels SET type=$1, snowflakes=$2, enabled=TRUE WHERE guild_id=$3 AND channel_id=$4;', log_type, sf if sf else None, ctx.guild.id, channel_obj.id)
                for e in current_entries:
                    if e['channel_id'] == channel_obj.id:
                        e.update({'type': log_type, 'snowflakes': sf if sf else None, 'enabled': True})
                msg = f'{self.emoji.get_random_emoji()} Log channel {channel_obj.mention} updated with type `{log_type or "general"}`.'
            else:
                await conn.execute('INSERT INTO statistic_channels (guild_id, channel_id, type, snowflakes, enabled) VALUES ($1, $2, $3, $4, TRUE);', ctx.guild.id, channel_obj.id, log_type, sf if sf else None)
                current_entries.append({
                    'guild_id': ctx.guild.id,
                    'channel_id': channel_obj.id,
                    'type': log_type,
                    'snowflakes': sf if sf else None,
                    'enabled': True
                })
                msg = f'{self.emoji.get_random_emoji()} Log channel {channel_obj.mention} created with type `{log_type or "general"}`.'
        await Statistics.load_channels()
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())

    # DONE
    @app_commands.command(name='rmv', description='Move all the members in one room to another.')
    @app_commands.describe(source_id='Tag the source channel or include its snowflake ID',target_id='Tag the target channel or include its snowflake ID')
    @is_owner_developer_administrator_predicator()
    async def room_move_all_app_command(
        self,
        interaction: discord.Interaction,
        source_id: Optional[str] = None,
        target_id: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        source_channel = interaction.guild.get_channel(int(source_id))
        target_channel = interaction.guild.get_channel(int(target_id))
        if not source_channel or not target_channel:
            return await interaction.response.send_message(content='\U0001F6AB One or both channels are invalid.')
        if not isinstance(source_channel, discord.VoiceChannel) or not isinstance(target_channel, discord.VoiceChannel):
            raise ValueError('\U0001F6AB Both source and target must be voice channels.')
        for member in source_channel.members:
            try:
                await member.move_to(target_channel)
            except discord.Forbidden:
                logger.warning(f'\U0001F6AB Missing permissions to move {member}.')
            except discord.HTTPException:
                logger.warning(f'\U0001F6AB Failed to move {member} due to a network error.')
        await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Moved all members from {source_channel.mention} to {target_channel.mention}.')
        
    # DONE
    @commands.command(name='rmv', help='Move all the members in one room to another.')
    @is_owner_developer_administrator_predicator()
    async def room_move_all_text_command(
        self,
        ctx: commands.Context,
        source_id: Optional[int] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        target_id: Optional[int] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        source_channel = ctx.guild.get_channel(source_id)
        target_channel = ctx.guild.get_channel(target_id)
        if not source_channel or not target_channel:
            await self.handler.send_message(ctx, content='\U0001F6AB One or both channel IDs are invalid.')
            return
        if not isinstance(source_channel, discord.VoiceChannel) or not isinstance(target_channel, discord.VoiceChannel):
            raise ValueError('\U0001F6AB Both source and target must be voice channels.')
        for member in source_channel.members:
            try:
                await member.move_to(target_channel)
            except discord.Forbidden:
                logger.warning(f'\U0001F6AB Missing permissions to move {member}.')
            except discord.HTTPException:
                logger.warning(f'\U0001F6AB Failed to move {member} due to a network error.')
        await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Moved all members from {source_channel.mention} to {target_channel.mention}.')
            
            
    # DONE
    @app_commands.command(name='smute', description='Mutes or unmutes a member throughout the entire guild.')
    @app_commands.describe(member='Tag a member or include their snowflake ID', reason='Optional reason (required for 7 days or more)')
    @is_owner_developer_administrator_predicator()
    async def toggle_server_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None,
        reason: Optional[str] = 'No reason provided'
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        member_obj = await self.member_service.resolve_member(interaction, member)
        if not member_obj:
            return await interaction.response.send_message(content=f'\U0001F6AB Could not resolve a valid member from target: `{member}`.')
        if member_obj.bot:
            return await interaction.response.send_message(content='\U0001F6AB You cannot server mute the bot.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT server_mute_guild_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            current_mutes = row.get('server_mute_guild_ids', []) if row else []
            action = None
            if not current_mutes:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, server_mute_guild_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET server_mute_guild_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(COALESCE (u.server_mute_guild_ids, '{}') || ARRAY[$2])
                        )
                        FROM users u WHERE u.discord_snowflake = EXCLUDED.discord_snowflake
                    ),
                    updated_at = NOW()
                ''', member_obj.id, interaction.guild.id)
                await conn.execute('''
                    INSERT INTO active_server_voice_mutes (guild_id, discord_snowflake, reason)
                    VALUES ($1, $2, $3)
                    ON CONFLICT (guild_id, discord_snowflake) DO UPDATE
                    SET reason = EXCLUDED.reason
                ''', interaction.guild.id, member_obj.id, reason)
                if member_obj.voice and member_obj.voice.channel:
                    try:
                        await member_obj.edit(mute=True)
                    except discord.Forbidden:
                        return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} was not successfully voice muted.', allowed_mentions=discord.AllowedMentions.none())
                action = 'muted'
            if interaction.guild.id in current_mutes:
                await conn.execute('UPDATE users SET server_mute_guild_ids = array_remove(server_mute_guild_ids, $2), updated_at = NOW() WHERE discord_snowflake = $1', member_obj.id, interaction.guild.id)
                await conn.execute('DELETE FROM active_server_voice_mutes WHERE guild_id = $1 AND discord_snowflake = $2', interaction.guild.id, member_obj.id)
                if member_obj.voice and member_obj.voice.channel:
                    try:
                        await member_obj.edit(mute=False)
                    except discord.Forbidden:
                        return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} was not successfully voice unmuted.', allowed_mentions=discord.AllowedMentions.none())
                action = 'unmuted'
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_server_mute', member_obj.id, interaction.user.id, interaction.guild.id, None, f'Server {action}')
        return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been server {action} for reason: {reason}', allowed_mentions=discord.AllowedMentions.none())
            
    # DONE
    @commands.command(name='smute', help='Mutes or unmutes a member throughout the entire guild.')
    @is_owner_developer_administrator_predicator()
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        *,
        reason: Optional[str] = commands.parameter(default='No reason provided', description='Optional reason (required for 7 days or more)')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot server mute the bot.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('SELECT server_mute_guild_ids FROM users WHERE discord_snowflake = $1', member_obj.id)
            current_mutes = row.get('server_mute_guild_ids', []) if row else []
            action = None
            if not current_mutes:
                await conn.execute('''
                    INSERT INTO users (discord_snowflake, server_mute_guild_ids)
                    VALUES ($1, ARRAY[$2]::BIGINT[])
                    ON CONFLICT (discord_snowflake) DO UPDATE
                    SET server_mute_guild_ids = (
                        SELECT ARRAY(
                            SELECT DISTINCT unnest(COALESCE (u.server_mute_guild_ids, '{}') || ARRAY[$2])
                        )
                        FROM users u WHERE u.discord_snowflake = EXCLUDED.discord_snowflake
                    ),
                    updated_at = NOW()
                ''', member_obj.id, ctx.guild.id)
                await conn.execute('''
                    INSERT INTO active_server_voice_mutes (guild_id, discord_snowflake, reason)
                    VALUES ($1, $2, $3)
                    ON CONFLICT (guild_id, discord_snowflake) DO UPDATE
                    SET reason = EXCLUDED.reason
                ''', ctx.guild.id, member_obj.id, reason)
                if member_obj.voice and member_obj.voice.channel:
                    try:
                        await member_obj.edit(mute=True)
                    except discord.Forbidden:
                        return await self.handler.send_message(content=f'\U0001F6AB {member_obj.mention} was not successfully voice muted.', allowed_mentions=discord.AllowedMentions.none())
                action = 'muted'
            if ctx.guild.id in current_mutes:
                await conn.execute('UPDATE users SET server_mute_guild_ids = array_remove(server_mute_guild_ids, $2), updated_at = NOW() WHERE discord_snowflake = $1', member_obj.id, ctx.guild.id)
                await conn.execute('DELETE FROM active_server_voice_mutes WHERE guild_id = $1 AND discord_snowflake = $2', ctx.guild.id, member_obj.id)
                if member_obj.voice and member_obj.voice.channel:
                    try:
                        await member_obj.edit(mute=False)
                    except discord.Forbidden:
                        return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} was not successfully voice unmuted.', allowed_mentions=discord.AllowedMentions.none())
                action = 'unmuted'
                await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_server_mute', member_obj.id, ctx.author.id, ctx.guild.id, None, f'Server {action}')
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been server {action} for reason: {reason}', allowed_mentions=discord.AllowedMentions.none())

    # DONE
    @app_commands.command(name='temps', description='List temporary rooms with matching command aliases.')
    @is_owner_developer_administrator_predicator()
    async def check_temp_rooms_app_command(self, interaction: discord.Interaction):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(interaction.guild)
        if not rooms:
            return await interaction.response.send_message(content='\U0001F6AB No temporary rooms found.')
        lines = []
        for room in rooms:
            room_id = room.room_snowflake
            channel_obj = await self.channel_service.resolve_channel(interaction, room_id)
            aliases = await Alias.fetch_command_aliases_by_channel_id(interaction.guild.id, room.room_snowflake)
            lines.append(f"{channel_obj.mention} ({room_id})")
            if not aliases:
                continue
            for alias_obj in aliases:
                lines.append(f"  ↳ {alias_obj.alias_name} ({alias_obj.alias_type})")
        pages = []
        chunks = 18
        for i in range(0, len(lines), chunks):
            step = lines[i:i + chunks]
            embed = discord.Embed(
                title=f"Temporary Rooms in {interaction.guild.name}",
                description='\n'.join(step),
                color=discord.Color.blue()
            )
            pages.append(embed)
        paginator = AppPaginator(self.bot, interaction, pages)
        await paginator.start()
        
    # DONE
    @commands.command(name='temps', help='List temporary rooms with matching command aliases.')
    @is_owner_developer_administrator_predicator()
    async def check_temp_rooms_text_command(self, ctx: commands.Context):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        rooms = await TemporaryRoom.fetch_temporary_rooms_by_guild(ctx.guild)
        if not rooms:
            return await self.handler.send_message(ctx, content='\U0001F6AB No temporary rooms found.')
        lines = []
        for room in rooms:
            room_id = room.room_snowflake
            channel_obj = await self.channel_service.resolve_channel(ctx, room_id)
            aliases = await Alias.fetch_command_aliases_by_channel_id(ctx.guild.id, room.room_snowflake)
            lines.append(f"{channel_obj.mention} ({room_id})")
            if not aliases:
                continue
            for alias_obj in aliases:
                lines.append(f"  ↳ {alias_obj.alias_name} ({alias_obj.alias_type})")
        pages = []
        chunk_size = 18
        for i in range(0, len(lines), chunk_size):
            step = lines[i:i + chunk_size]
            embed = discord.Embed(
                title=f"Temporary Rooms in {ctx.guild.name}",
                description='\n'.join(step),
                color=discord.Color.blue()
            )
            pages.append(embed)
        paginator = Paginator(self.bot, ctx, pages)
        await paginator.start()
        
    # DONE
    @app_commands.command(name='xcap', description='Delete a duration cap for bans, mutes and text mutes.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(channel='Tag a channel or include its snowflake ID', moderation_type='One of: `mute`, `ban`, `tmute`')
    async def undo_cap_interaction(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None,
        moderation_type: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        valid_types = {'mute','ban','tmute'}
        if moderation_type not in valid_types:
            return await interaction.response.send_message(content=f'\U0001F6AB Invalid moderation type. Must be one of: {", ".join(valid_types)}')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4',
                interaction.guild.id, channel_obj.id, moderation_type, channel_obj.name
            )
            if row:
                original_duration = row['duration']
                msg = f'{self.emoji.get_random_emoji()} Cap deleted on {channel_obj.mention} for {moderation_type} of duration {original_duration}.'
            else:
                return await interaction.response.send_message(content=f'\U0001F6AB No caps exist in {channel_obj.mention} for {moderation_type}.')
            await conn.execute('''
                DELETE FROM active_caps
                WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4
            ''', interaction.guild.id, channel_obj.id, moderation_type, channel_obj.name)
        await interaction.response.send_message(content=msg, allowed_mentions=discord.AllowedMentions.none())
    
    # DONE
    @commands.command(name='xcap', help='Delete a duration cap for bans, mutes and text mutes.')
    @is_owner_developer_administrator_predicator()
    async def undo_cap_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        moderation_type: Optional[str] = commands.parameter(default=None, description='One of: `mute`, `ban`, `tmute`')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        valid_types = {'mute', 'ban', 'tmute'}
        if moderation_type not in valid_types:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Invalid moderation type. Must be one of: {", ".join(valid_types)}')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow(
                'SELECT duration FROM active_caps WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4',
                ctx.guild.id, channel_obj.id, moderation_type, channel_obj.name
            )
            if row:
                original_duration = row['duration']
                msg = f'{self.emoji.get_random_emoji()} Cap deleted on {channel_obj.mention} for {moderation_type} of duration {original_duration}.'
            else:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No caps exist in {channel_obj.mention} for {moderation_type}.')
            await conn.execute('''
                DELETE FROM active_caps
                WHERE guild_id=$1 AND channel_id=$2 AND moderation_type=$3 AND room_name=$4
            ''', ctx.guild.id, channel_obj.id, moderation_type, channel_obj.name)
        await self.handler.send_message(ctx, content=msg, allowed_mentions=discord.AllowedMentions.none())
        
    # DONE
    @app_commands.command(name='xstage', description='Destroy the stage in the current channel.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_predicator()
    async def stage_quit_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content=f'\U0001F6AB {channel_obj.mention} is not a valid target.')
        async with self.bot.db_pool.acquire() as conn:
            stage = await conn.fetchrow('SELECT initiator_id FROM active_stages WHERE guild_id=$1 AND channel_id=$2', interaction.guild.id, channel_obj.id)
            if not stage:
                return await interaction.response.send_message(content='\U000026A0\U0000FE0F No active stage found.')
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content='\U0001F6AB Only the admins and above can end this stage.')
            await conn.execute('''
                DELETE FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', interaction.guild.id, channel_obj.id, channel_obj.name)
            await conn.execute('''
                DELETE FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', interaction.guild.id, channel_obj.id, channel_obj.name)
            await conn.execute('DELETE FROM active_voice_mutes WHERE guild_id=$1 AND channel_id=$2 AND target=\'room\'', interaction.guild.id, channel_obj.id)
            for member in channel_obj.members:
                user_mute = await conn.fetchrow('''
                    SELECT 1 FROM active_voice_mutes
                    WHERE guild_id=$1 AND channel_id=$2 AND discord_snowflake=$3 AND target='room' AND room_name=$4
                ''', interaction.guild.id, channel_obj.id, member.id, channel_obj.name)
                if not user_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(mute=False, reason='Stage ended — no user-specific mute found')
                    except discord.Forbidden:
                        logger.warning(f'Failed to mute member: {member.mention}.')
        return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Stage in {channel_obj.mention} has ended.')
        
    # DONE
    @commands.command(name='xstage', help='Destroy the stage in the current channel.')
    @is_owner_developer_administrator_predicator()
    async def stage_quit_text_command(
        self,
        ctx: commands.Context,
        *,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB {channel_obj.mention} is not a valid target.')
        async with self.bot.db_pool.acquire() as conn:
            stage = await conn.fetchrow('SELECT initiator_id FROM active_stages WHERE guild_id=$1 AND channel_id=$2', ctx.guild.id, channel_obj.id)
            if not stage:
                return await self.handler.send_message(ctx, content='\U000026A0\U0000FE0F No active stage found.')
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB Only the admins and above can end this stage.')
            await conn.execute('''
                DELETE FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', ctx.guild.id, channel_obj.id, channel_obj.name)
            await conn.execute('''
                DELETE FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', ctx.guild.id, channel_obj.id, channel_obj.name)
            await conn.execute('DELETE FROM active_voice_mutes WHERE guild_id=$1 AND channel_id=$2 AND target=\'room\'', ctx.guild.id, channel_obj.id)
            for member in channel_obj.members:
                user_mute = await conn.fetchrow('''
                    SELECT 1 FROM active_voice_mutes
                    WHERE guild_id=$1 AND channel_id=$2 AND discord_snowflake=$3 AND target='room' AND room_name=$4
                ''', ctx.guild.id, channel_obj.id, member.id, channel_obj.name)
                if not user_mute and member.voice and member.voice.mute:
                    try:
                        await member.edit(mute=False, reason='Stage ended — no user-specific mute found')
                    except discord.Forbidden:
                        logger.warning(f'Failed to unmute {member.mention}')
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Stage in {channel_obj.mention} has ended.')

async def setup(bot: DiscordBot):
    await bot.add_cog(AdminCommands(bot))
