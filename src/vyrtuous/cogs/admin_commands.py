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
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.check_service import *
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.role_service import RoleService
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.coordinator import Coordinator
from vyrtuous.utils.duration import AppDuration, Duration
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.moderator import Moderator
from vyrtuous.utils.server_mute import ServerMute
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.state import State
from vyrtuous.utils.statistics import Statistics
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.voice_mute import VoiceMute
from vyrtuous.service.message_service import MessageService

class AdminCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.handler = MessageService(self.bot, self.bot.db_pool)
        self.member_service = MemberService()
        self.role_service = RoleService()
    
    # DONE
    @app_commands.command(name='alias', description='Set an alias for a Vyrtuous action.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(
        alias_name='Alias/Pseudonym',
        alias_type='One of: cow, uncow, mute, unmute, ban, unban, flag, unflag, tmute, untmute, role, unrole',
        channel='Tag a channel or include its snowflake ID',
        role='Role ID (only for role/unrole)'
    )
    async def create_alias_app_command(
        self,
        interaction: discord.Interaction,
        alias_name: Optional[str],
        alias_type: Optional[str],
        channel: AppChannelSnowflake,
        role: AppRoleSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        try:
            role_obj = await self.role_service.resolve_role(interaction, role)
            role_snowflake = role_obj.id
        except Exception as e:
            role_snowflake = None
        alias = Alias(alias_name=alias_name, alias_type=alias_type, channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, role_snowflake=role_snowflake)
        await alias.grant()
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_alias', None, interaction.user.id, interaction.guild.id, channel_obj.id, f'Created an alias: {alias_name}')
        if role:
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Alias `{alias_name}` created successfully for channel {channel_obj.mention} and role {role_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {str(e)}')
        else:
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Alias `{alias_name}` created successfully for channel {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @commands.command(name='alias', help='Set an alias for a cow, uncow, mute, unmute, ban, unban, flag, unflag, tmute, untmute, role, or unrole action.')
    @is_owner_developer_administrator_predicator()
    async def create_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: Optional[str] = commands.parameter(default=None, description='Alias/Pseudonym'),
        alias_type: Optional[str] = commands.parameter(default=None, description='One of: `cow`, `uncow`, `mute`, `unmute`, `ban`, `unban`, `flag`, `unflag`, `tmute`, `untmute`, `role`, `unrole`'),
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        *,
        role: RoleSnowflake = commands.parameter(default=None, description='Role ID (only for role/unrole)')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        try:
            role_obj = await self.role_service.resolve_role(ctx, role)
            role_snowflake = role_obj.id
        except Exception as e:
            role_snowflake = None
        alias = Alias(alias_name=alias_name, alias_type=alias_type, channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, role_snowflake=role_snowflake)
        await alias.grant()
        async with ctx.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'create_alias', None, ctx.author.id, ctx.guild.id, channel_obj.id, f'Created an alias: {alias_name}')
        if role:
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Alias `{alias_name}` created successfully for channel {channel_obj.mention} and role {role_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {str(e)}')
        else:
            try:
                return await state.end(success=f'{self.emoji.get_random_emoji()} Alias `{alias_name}` created successfully for channel {channel_obj.mention}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {str(e)}')
    
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
        duration: AppDuration = '24h'
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        duration_seconds = duration.build_timedelta().total_seconds()
        cap = await Cap.fetch_by_channel_guild_and_moderation_type(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, moderation_type=moderation_type)
        if cap:
            await Cap.update_by_channel_and_duration(channel_snowflake=channel_obj.id, duration=duration_seconds)
        else:
            cap = Cap(channel_snowflake=channel_obj.id, duration=duration_seconds, guild_snowflake=interaction.guild.id, moderation_type=moderation_type)
            await cap.grant()
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Cap `{moderation_type}` created for {channel_obj.mention} successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @commands.command(name='cap', help='Set a duration limit for bans, mutes and text mutes.')
    @is_owner_developer_administrator_predicator()
    async def cap_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        moderation_type: Optional[str] = commands.parameter(default=None, description='One of: `mute`, `ban`, `tmute`'),
        *,
        duration: Duration = commands.parameter(default='24h', description='Options: (+|-)duration(m|h|d) \n 0 - permanent / 24h - default')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        duration_seconds = duration.build_timedelta().total_seconds()
        cap = await Cap.fetch_by_channel_guild_and_moderation_type(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, moderation_type=moderation_type)
        if cap:
            await Cap.update_by_channel_and_duration(channel_snowflake=channel_obj.id, duration=duration_seconds)
        else:
            cap = Cap(channel_snowflake=channel_obj.id, duration=duration_seconds, guild_snowflake=ctx.guild.id, moderation_type=moderation_type)
            await cap.grant()
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Cap `{moderation_type}` created for {channel_obj.mention} successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @app_commands.command(name='chown', description='Change the owner of a temporary room.')
    @app_commands.describe(member='Tag a user or provide their snowflake ID', channel='Tag a channel or provide it\'s snowflake ID')
    @is_owner_developer_administrator_predicator()
    async def change_temp_room_owner_app_command(
        self,
        interaction,
        channel: AppChannelSnowflake,
        member: AppMemberSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
            member_obj = await self.member_service.resolve_member(interaction, member)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        await TemporaryRoom.update_owner(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room `{channel_obj.mention} ownership transfered to {member_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
            
    # DONE
    @commands.command(name='chown', help='Change the owner of a temporary room.')
    @is_owner_developer_administrator_predicator()
    async def change_temp_room_owner_text_command(
        self,
        ctx,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or provide it\'s snowflake ID'),
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a user or provide their snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
            member_obj = await self.member_service.resolve_member(ctx, member)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        await TemporaryRoom.update_owner(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room `{channel_obj.mention} ownership transfered to {member_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @app_commands.command(name='coord', description='Grants/revokes coordinator access for a specific voice channel.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(member='Tag a member or include their snowflake ID', channel='Tag a channel or include its snowflake ID')
    async def create_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
            await has_equal_or_higher_role(interaction, channel_snowflake=channel_obj.id, member_snowflake=member_obj.id)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        channel_snowflakes = await Coordinator.fetch_channels_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
        if channel_snowflakes and channel_obj.id in channel_snowflakes:
            await Coordinator.delete_by_channel_and_member(channel_snowflake=channel_obj.id, member_snowflake=member_obj.id)
            action = 'revoked'
        else:
            coordinator = Coordinator(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            await coordinator.grant()
            action = 'granted'
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_coordinator', member_obj.id, interaction.user.id, interaction.guild.id, channel_obj.id, f'Coordinator access {action}')
        try:
            return await state.end(success=f"{self.emoji.get_random_emoji()} Cooridinator access has been {action} for {member_obj.mention} in {channel_obj.mention}.")
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @commands.command(name='coord', help='Grants/revokes coordinator access for a specific voice channel.')
    @is_owner_developer_administrator_predicator()
    async def create_coordinator_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
            await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
            await has_equal_or_higher_role(ctx, channel_snowflake=channel_obj.id, member_snowflake=member_obj.id)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        channel_snowflakes = await Coordinator.fetch_channels_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
        if channel_snowflakes and channel_obj.id in channel_snowflakes:
            await Coordinator.delete_by_channel_and_member(channel_snowflake=channel_obj.id, member_snowflake=member_obj.id)
            action = 'revoked'
        else:
            coordinator = Coordinator(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            await coordinator.grant()
            action = 'granted'
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_coordinator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, f'Coordinator access {action}')
        try:
            return await state.end(success=f"{self.emoji.get_random_emoji()} Cooridinator access has been {action} for {member_obj.mention} in {channel_obj.mention}.")
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
   # DONE
    @app_commands.command(name='cstage', description='Create a stage in the current or specified channel.')
    @app_commands.describe(channel='Tag a voice/stage channel', duration='Duration of the stage (e.g., 1h, 30m)')
    @is_owner_developer_administrator_predicator()
    async def stage_create_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        duration: AppDuration = 1
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        expires_at = duration.output_datetime()
        duration_display = duration.output_display()
        skipped, muted, failed = [], [], []
        stage = Stage(channel_snowflake=channel_obj.id, expires_at=expires_at, guild_snowflake=interaction.guild.id, member_snowflake=interaction.user.id)
        await stage.grant()
        for member in channel_obj.members:
            if is_owner_developer_administrator_coordinator_moderator_via_channel_member(channel_snowflake=channel_obj.id, member_snowflake=member.id) or member.id == interaction.user.id:
                skipped.append(member)
                continue
            voice_mute = await VoiceMute(channel_snowflake=channel_obj.id, expires_at=expires_at, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target="room", reason="Stage mute")
            await voice_mute.grant()
            try:
                if member.voice and member.voice.channel.id == channel_obj.id:
                    await member.edit(mute=True)
                muted.append(member)
            except Exception as e:
                logger.warning(f'Failed to mute.')
                failed.append(member)
        description_lines = [
            f"**Duration:** {duration_display}",
            f"**Channel:** {channel_obj.mention}",
            f"**Muted:** {len(muted)} user(s)",
            f"**Skipped:** {len(skipped)}"
        ]
        if failed:
            description_lines.append(f"**Failed:** {len(failed)}")
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{self.emoji.get_random_emoji()} Stage Created",
            color=discord.Color.blurple()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @commands.command(name='cstage', help='Create a stage in the current or specified channel.')
    @is_owner_developer_administrator_predicator()
    async def stage_create_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        *,
        duration: Duration = '1'
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
            await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        expires_at = duration.output_datetime()
        duration_display = duration.output_display()
        skipped, muted, failed = [], [], []
        stage = Stage(channel_snowflake=channel_obj.id, expires_at=expires_at, guild_snowflake=ctx.guild.id, member_snowflake=ctx.author.id)
        await stage.grant()
        for member in channel_obj.members:
            if is_owner_developer_administrator_coordinator_moderator_via_channel_member(channel_snowflake=channel_obj.id, member_snowflake=member.id) or member.id == ctx.author.id:
                skipped.append(member)
                continue
            voice_mute = await VoiceMute(channel_snowflake=channel_obj.id, expires_at=expires_at, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target="room", reason="Stage mute")
            await voice_mute.grant()
            try:
                if member.voice and member.voice.channel.id == channel_obj.id:
                    await member.edit(mute=True)
                muted.append(member)
            except Exception as e:
                logger.warning(f'Failed to mute.')
                failed.append(member)
        description_lines = [
            f"**Duration:** {duration_display}",
            f"**Channel:** {channel_obj.mention}",
            f"**Muted:** {len(muted)} user(s)",
            f"**Skipped:** {len(skipped)}"
        ]
        if failed:
            description_lines.append(f"**Failed:** {len(failed)}")
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{self.emoji.get_random_emoji()} Stage Created",
            color=discord.Color.blurple()
        )
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @app_commands.command(name='log', description='Toggle logging for a channel on or off.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    async def toggle_log_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        statistics = Statistics.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        if statistics:
            enabled = "start"
        else:
            enabled = "stop"
        await Statistics.update_by_channel_enabled_and_guild(channel_snowflake=channel_obj.id, enabled=enabled, guild_snowflake=interaction.guild.id)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Logging messages will now {enabled} being sent to {channel_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @commands.command(name='log', help='Toggle logging for a channel on or off.')
    @is_owner_developer_administrator_predicator()
    async def toggle_log_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
            await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        statistics = Statistics.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        if statistics:
            enabled = "start"
        else:
            enabled = "stop"
        await Statistics.update_by_channel_enabled_and_guild(channel_snowflake=channel_obj.id, enabled=enabled, guild_snowflake=ctx.guild.id)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Logging messages will now {enabled} being sent to {channel_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @app_commands.command(name='logs', description='Lists log channels for this guild or a specified guild ID.')
    @is_owner_developer_administrator_predicator()
    async def list_logs_app_command(
        self,
        interaction: discord.Interaction
    ):
        state = State(interaction)
        statistics = await Statistics.fetch_by_guild(guild_snowflake=interaction.guild.id)
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Log Channels',
            color=discord.Color.blue()
        )
        embed.set_footer(text=f'Guild ID: {interaction.guild.id}')
        for statistic in statistics:
            enabled = statistic.enabled
            snowflakes = statistic.snowflakes
            statistic_type = statistic.statistic_type
            if statistic_type == 'general':
                detail = f"Logs all events in {interaction.guild.name}"
            elif statistic_type == 'channel':
                detail=f"Logs only events from channels: {', '.join(f'<#{s}>' for s in snowflakes)}"
            elif statistic_type == 'member':
                detail=f"Logs only events for members: {', '.join(f'<@{s}>' for s in snowflakes)}"
            embed.add_field(name=f"{'✅' if enabled else '\U000026A0\U0000FE0F'}", value=f"Type: **{statistic_type}**\n{detail}", inline=False)
        try:
            return await state.end(embed=embed) 
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
               
    # DONE
    @commands.command(name='logs', help='Lists log channels for this guild or a specified guild ID.')
    @is_owner_developer_administrator_predicator()
    async def list_logs_text_command(
        self,
        ctx: commands.Context
    ):
        state = State(ctx)
        statistics = await Statistics.fetch_by_guild(guild_snowflake=ctx.guild.id)
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Log Channels',
            color=discord.Color.blue()
        )
        embed.set_footer(text=f'Guild ID: {ctx.guild.id}')
        for statistic in statistics:
            enabled = statistic.enabled
            snowflakes = statistic.snowflakes
            statistic_type = statistic.statistic_type
            if statistic_type == 'general':
                detail = f"Logs all events in {ctx.guild.name}"
            elif statistic_type == 'channel':
                detail=f"Logs only events from channels: {', '.join(f'<#{s}>' for s in snowflakes)}"
            elif statistic_type == 'member':
                detail=f"Logs only events for members: {', '.join(f'<@{s}>' for s in snowflakes)}"
            embed.add_field(name=f"{'✅' if enabled else '\U000026A0\U0000FE0F'}", value=f"Type: **{statistic_type}**\n{detail}", inline=False)
        try:
            return await state.end(embed=embed) 
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')

    # DONE
    @app_commands.command(name='mlog', description='Create, modify, or delete a log channel.')
    @app_commands.describe(
        channel='Tag a channel or include its snowflake ID',
        action='create | modify | delete',
        statistic_type='Type of logs: member, channel, general',
        snowflakes='Optional list of member IDs to include in logs'
    )
    @is_owner_developer_administrator_predicator()
    async def modify_log_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        action: Optional[str] = None,
        statistic_type: Optional[str] = None,
        snowflakes: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel) 
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        match action.lower():
            case 'create':
                statistics = Statistics(channel_snowflake=channel_obj.id, enabled=True, guild_snowflake=interaction.guild.id, statistic_type=statistic_type, snowflakes=snowflakes)
                await statistics.grant()
                action = "created"
            case 'delete':
                await Statistics.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
                action = "deleted"
            case 'modify':
                await Statistics.update_by_channel_guild_and_type(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, statistic_type=statistic_type, snowflakes=snowflakes)
                action = "modified"
        try:
            return await state.end(f'{self.emoji.get_random_emoji()} Log {action} successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
        
    # DONE
    @commands.command(name='mlog', help='Create, modify, or delete a log channel.')
    @is_owner_developer_administrator_predicator()
    async def modify_log_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        action: Optional[str] = commands.parameter(default=None, description='create | modify | delete'),
        statistic_type: Optional[str] = commands.parameter(default=None, description='Type of logs: member, channel, general'),
        *snowflakes: Optional[int]
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel) 
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        match action.lower():
            case 'create':
                statistics = Statistics(channel_snowflake=channel_obj.id, enabled=True, guild_snowflake=ctx.guild.id, statistic_type=statistic_type, snowflakes=snowflakes)
                await statistics.grant()
                action = "created"
            case 'delete':
                await Statistics.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
                action = "deleted"
            case 'modify':
                await Statistics.update_by_channel_guild_and_type(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, statistic_type=statistic_type, snowflakes=snowflakes)
                action = "modified"
        try:
            return await state.end(f'{self.emoji.get_random_emoji()} Log {action} successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')

    # DONE
    @app_commands.command(name='rmv', description='Move all the members in one room to another.')
    @app_commands.describe(
        source_channel='Tag the source channel or include its snowflake ID',
        target_channel='Tag the target channel or include its snowflake ID'
    )
    @is_owner_developer_administrator_predicator()
    async def room_move_all_app_command(
        self,
        interaction: discord.Interaction,
        source_channel: AppChannelSnowflake,
        target_channel: AppChannelSnowflake
    ):
        state = State(interaction)
        moved = []
        failed = []
        try:
            source_channel_obj = await self.channel_service.resolve_channel(interaction, source_channel)
            target_channel_obj = await self.channel_service.resolve_channel(interaction, target_channel)
        except:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        for member in source_channel_obj.members:
            try:
                await member.move_to(target_channel_obj)
                moved.append(member)
            except discord.Forbidden:
                failed.append(member)
                logger.warning(f'\U000026A0\U0000FE0F Missing permissions to move {member}.')
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Room Move Results',
            color=discord.Color.green()
        )
        if moved:
            embed.add_field(name=f'Successfully Moved ({len(moved)})', value='\n'.join(member.mention for member in moved), inline=False)
        else:
            embed.add_field(name='Successfully Moved', value='None', inline=False)
        if failed:
            embed.add_field(name=f'Failed to Move ({len(failed)})', value='\n'.join(member.mention for member in failed), inline=False)
        embed.set_footer(text=f'Moved from {source_channel_obj.name} to {target_channel_obj.name}')
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
        
    # DONE
    @commands.command(name='rmv', help='Move all the members in one room to another.')
    @is_owner_developer_administrator_predicator()
    async def room_move_all_text_command(
        self,
        ctx: commands.Context,
        source_channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        target_channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        moved = []
        failed = []
        try:
            source_channel_obj = await self.channel_service.resolve_channel(ctx, source_channel)
            target_channel_obj = await self.channel_service.resolve_channel(ctx, target_channel)
        except:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        for member in source_channel_obj.members:
            try:
                await member.move_to(target_channel_obj)
                moved.append(member)
            except discord.Forbidden:
                failed.append(member)
                logger.warning(f'\U000026A0\U0000FE0F Missing permissions to move {member}.')
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Room Move Results',
            color=discord.Color.green()
        )
        if moved:
            embed.add_field(name=f'Successfully Moved ({len(moved)})', value='\n'.join(member.mention for member in moved), inline=False)
        else:
            embed.add_field(name='Successfully Moved', value='None', inline=False)
        if failed:
            embed.add_field(name=f'Failed to Move ({len(failed)})', value='\n'.join(member.mention for member in failed), inline=False)
        embed.set_footer(text=f'Moved from {source_channel_obj.name} to {target_channel_obj.name}')
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
            
    # DONE
    @app_commands.command(name='smute', description='Mutes or unmutes a member throughout the entire guild.')
    @app_commands.describe(member='Tag a member or include their snowflake ID', reason='Optional reason (required for 7 days or more)')
    @is_owner_developer_administrator_predicator()
    async def toggle_server_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        reason: Optional[str] = 'No reason provided'
    ):
        state = State(interaction)
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(interaction, member)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        if member_obj.id == interaction.guild.me.id:
            return await state.end(warning="You cannot server mute {interaction.guild.me.mention}")
        server_mute = await ServerMute.fetch_by_member(member_snowflake=member_obj.id)
        if not server_mute:
            server_mute = ServerMute(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, reason=reason)
            server_mute.grant()
            action = 'muted'
            should_be_muted = True
        else:
            await ServerMute.delete_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            action = 'unmuted'
            should_be_muted = False
        if member_obj.voice and member_obj.voice.channel:
            try:
                await member_obj.edit(mute=should_be_muted)
            except discord.Forbidden as e:
                return await state.end(error=f'\U0001F3C6 {str(e)}')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_server_mute', member_obj.id, interaction.user.id, interaction.guild.id, None, f'Server {action}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully server {action} {member_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
            
    # DONE
    @commands.command(name='smute', help='Mutes or unmutes a member throughout the entire guild.')
    @is_owner_developer_administrator_predicator()
    async def toggle_server_mute_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        *,
        reason: Optional[str] = commands.parameter(default='No reason provided', description='Optional reason (required for 7 days or more)')
    ):
        state = State(ctx)
        member_obj = None
        try:
            member_obj = await self.member_service.resolve_member(ctx, member)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        if member_obj.id == ctx.guild.me.id:
            return await state.end(warning="You cannot server mute {ctx.guild.me.mention}")
        server_mute = await ServerMute.fetch_by_member(member_snowflake=member_obj.id)
        if not server_mute:
            server_mute = ServerMute(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, reason=reason)
            server_mute.grant()
            action = 'muted'
            should_be_muted = True
        else:
            await ServerMute.delete_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            action = 'unmuted'
            should_be_muted = False
        if member_obj.voice and member_obj.voice.channel:
            try:
                await member_obj.edit(mute=should_be_muted)
            except discord.Forbidden as e:
                return await state.end(error=f'\U0001F3C6 {str(e)}')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_server_mute', member_obj.id, ctx.author.id, ctx.guild.id, None, f'Server {action}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Successfully server {action} {member_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
        
    # DONE
    @app_commands.command(name='temp', description='Toggle a temporary room and assign an owner.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID', owner='Tag a member or include their snowflake ID')
    @is_owner_developer_administrator_predicator()
    async def toggle_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        owner: AppMemberSnowflake
    ):
        state = State(interaction)
        action = None
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        temporary_room = await TemporaryRoom.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        if temporary_room:
            if temporary_room.member_snowflake:
                await Moderator.delete_by_channel_and_member(channel_snowflake=channel_obj.id, member_snowflake=temporary_room.member_snowflake)
            await TemporaryRoom.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            await Alias.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            action = 'removed'
        else:
            try:
                member_obj = await self.member_service.resolve_member(interaction, owner)
            except Exception as e:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
            temporary_room = TemporaryRoom(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, room_name=channel_obj.name)
            await temporary_room.grant()
            moderator = Moderator(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            await moderator.grant()
            action = f'created'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room {action}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
        
    # DONE
    @commands.command(name='temp', help='Toggle a temporary room and assign an owner.')
    @is_owner_developer_administrator_predicator()
    async def toggle_temp_room_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        owner: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their Discord ID')
    ):
        state = State(ctx)
        action = None
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        temporary_room = await TemporaryRoom.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        if temporary_room:
            if temporary_room.member_snowflake:
                await Moderator.delete_by_channel_and_member(channel_snowflake=channel_obj.id, member_snowflake=temporary_room.member_snowflake)
            await TemporaryRoom.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            await Alias.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            action = 'removed'
        else:
            try:
                member_obj = await self.member_service.resolve_member(ctx, owner)
            except Exception as e:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
            temporary_room = TemporaryRoom(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, room_name=channel_obj.name)
            await temporary_room.grant()
            moderator = Moderator(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            await moderator.grant()
            action = f'created'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Temporary room {action}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
        
    # DONE
    @app_commands.command(name='temps', description='List temporary rooms with matching command aliases.')
    @is_owner_developer_administrator_predicator()
    async def check_temp_rooms_app_command(self, interaction: discord.Interaction):
        state = State(interaction)
        pages = []
        channel_obj = None
        rooms = await TemporaryRoom.fetch_by_guild(guild_snowflake=interaction.guild.id)
        if not rooms:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms found.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}')
        lines = []
        for room in rooms:
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=room.channel_snowflake, guild_snowflake=interaction.guild.id)
            try:
                channel_obj = await self.channel_service.resolve_channel(interaction, room.channel_snowflake)
            except Exception as e:
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
                except Exception as e:
                    return await state.end(error=f'\U0001F3C6 {e}')
            lines.append(f"{channel_obj.mention} ({room.channel_snowflake})")
            if aliases:
                for alias in aliases:
                    lines.append(f"  ↳ {alias.alias_name} ({alias.alias_type})")
        chunks = 18
        for i in range(0, len(lines), chunks):
            step = lines[i:i + chunks]
            embed = discord.Embed(
                title=f"{self.emoji.get_random_emoji()} Temporary Rooms",
                description='\n'.join(step),
                color=discord.Color.blue()
            )
            pages.append(embed)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
        
    # DONE
    @commands.command(name='temps', help='List temporary rooms with matching command aliases.')
    @is_owner_developer_administrator_predicator()
    async def check_temp_rooms_text_command(self, ctx: commands.Context):
        state = State(ctx)
        pages = []
        channel_obj = None
        rooms = await TemporaryRoom.fetch_by_guild(guild_snowflake=ctx.guild.id)
        if not rooms:
            return await state.end(warning=f'\U000026A0\U0000FE0F No temporary rooms found in {ctx.guild.name}.')
        lines = []
        for room in rooms:
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=room.channel_snowflake, guild_snowflake=ctx.guild.id)
            try:
                channel_obj = await self.channel_service.resolve_channel(ctx, room.channel_snowflake)
            except Exception as e:
                return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
            lines.append(f"{channel_obj.mention} ({room.channel_snowflake})")
            if aliases:
                for alias in aliases:
                    lines.append(f"  ↳ {alias.alias_name} ({alias.alias_type})")
        chunks = 18
        for i in range(0, len(lines), chunks):
            step = lines[i:i + chunks]
            embed = discord.Embed(
                title=f"{self.emoji.get_random_emoji()} Temporary Rooms",
                description='\n'.join(step),
                color=discord.Color.blue()
            )
            pages.append(embed)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
        
    # DONE
    @app_commands.command(name='xalias', description='Deletes an alias.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(alias_name='Include an alias name')
    async def delete_alias_app_command(
        self,
        interaction: discord.Interaction,
        alias_name: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        alias = await Alias.fetch_by_guild_and_name(alias_name=alias_name, guild_snowflake=interaction.guild.id)
        if not alias:
            return await state.end(warning=f'{self.emoji.get_random_emoji()} No aliases found for `{alias_name}`.')
        await Alias.delete_by_guild_and_name(alias_name=alias.alias_name, guild_snowflake=interaction.guild.id)
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, alias.channel_snowflake)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1,$2,$3,$4,$5,$6)
                ''', 'delete_alias', None, interaction.user.id, interaction.guild.id, channel_obj.id, f'Deleted alias {alias_name}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Alias `{alias.alias_name}` of type `{alias.alias_type}` for channel {channel_obj.mention} deleted successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')

    # DONE
    @commands.command(name='xalias', help='Deletes an alias.')
    @is_owner_developer_administrator_predicator()
    async def delete_alias_text_command(
        self,
        ctx: commands.Context,
        alias_name: Optional[str] = commands.parameter(default=None, description='Include an alias name')
    ):
        state = State(ctx)
        channel_obj = None
        alias = await Alias.fetch_by_guild_and_name(alias_name=alias_name, guild_snowflake=ctx.guild.id)
        if not alias:
            return await state.end(warning=f'{self.emoji.get_random_emoji()} No aliases found for `{alias_name}`.')
        await Alias.delete_by_guild_and_name(alias_name=alias.alias_name, guild_snowflake=ctx.guild.id)
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, alias.channel_snowflake)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1,$2,$3,$4,$5,$6)
                ''', 'delete_alias', None, ctx.author.id, ctx.guild.id, channel_obj.id, f'Deleted alias {alias_name}')
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Alias `{alias.alias_name}` of type `{alias.alias_type}` for channel {channel_obj.mention} deleted successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @app_commands.command(name='xcap', description='Delete a duration cap for bans, mutes and text mutes.')
    @is_owner_developer_administrator_predicator()
    @app_commands.describe(channel='Tag a channel or include its snowflake ID', moderation_type='One of: `mute`, `ban`, `tmute`')
    async def undo_cap_interaction(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake,
        moderation_type: Optional[str] = None
    ):
        state = State(interaction)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        await Cap.delete_by_channel_guild_moderation_type(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, moderation_type=moderation_type)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Cap of type {moderation_type} and channel {channel_obj.mention} deleted successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @commands.command(name='xcap', help='Delete a duration cap for bans, mutes and text mutes.')
    @is_owner_developer_administrator_predicator()
    async def undo_cap_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID'),
        moderation_type: Optional[str] = commands.parameter(default=None, description='One of: `mute`, `ban`, `tmute`')
    ):
        state = State(ctx)
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {str(e)}')
        await Cap.delete_by_channel_guild_moderation_type(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, moderation_type=moderation_type)
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Cap of type {moderation_type} and channel {channel_obj.mention} deleted successfully.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @app_commands.command(name='xstage', description='Destroy the stage in the current channel.')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_predicator()
    async def stage_quit_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        failed_members = []
        succeeded_members = []
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except:
            channel_obj = interaction.channel
            await self.handler.send_message(interaction, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        stage = await Stage.fetch_by_guild_and_channel(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        if not stage:
            await state.end(warning=f'{self.emoji.get_random_emoji()} No stage for {channel_obj.mention} found.')
        await Stage.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        await VoiceMute.delete_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target="room")
        for member in channel_obj.members:
            voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target="user")
            if not voice_mute and member.voice and member.voice.mute:
                try:
                    await member.edit(mute=False, reason='Stage ended — no user-specific mute found')
                    succeeded_members.append(member)
                except discord.Forbidden:
                    logger.warning(f'Failed to unmute member: {member}')
                    failed_members.append(member)
        embed = discord.Embed(
            title="Stage Ended",
            color=discord.Color.green()
        )
        if succeeded_members:
            embed.add_field(name=f"Successfully Unmuted ({len(succeeded_members)})", value='\n'.join(member.mention for member in succeeded_members), inline=False
            )
        if failed_members:
            embed.add_field(name=f"Failed to Unmute ({len(failed_members)})", value='\n'.join(member.mention for member in failed_members), inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
    # DONE
    @commands.command(name='xstage', help='Destroy the stage in the current channel.')
    @is_owner_developer_administrator_predicator()
    async def stage_quit_text_command(
        self,
        ctx: commands.Context,
        *,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        failed_members = []
        succeeded_members = []
        channel_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except:
            channel_obj = ctx.channel
            await self.handler.send_message(ctx, content=f'\U000026A0\U0000FE0F Defaulting to {channel_obj.mention}.')
        stage = await Stage.fetch_by_guild_and_channel(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        if not stage:
            await state.end(warning=f'{self.emoji.get_random_emoji()} No stage for {channel_obj.mention} found.')
        await Stage.delete_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        await VoiceMute.delete_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target="room")
        for member in channel_obj.members:
            voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target="user")
            if not voice_mute and member.voice and member.voice.mute:
                try:
                    await member.edit(mute=False, reason='Stage ended — no user-specific mute found')
                    succeeded_members.append(member)
                except discord.Forbidden:
                    logger.warning(f'Failed to unmute member: {member}')
                    failed_members.append(member)
        embed = discord.Embed(
            title="Stage Ended",
            color=discord.Color.green()
        )
        if succeeded_members:
            embed.add_field(name=f"Successfully Unmuted ({len(succeeded_members)})", value='\n'.join(member.mention for member in succeeded_members), inline=False
            )
        if failed_members:
            embed.add_field(name=f"Failed to Unmute ({len(failed_members)})", value='\n'.join(member.mention for member in failed_members), inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {str(e)}')
        
async def setup(bot: DiscordBot):
    await bot.add_cog(AdminCommands(bot))
