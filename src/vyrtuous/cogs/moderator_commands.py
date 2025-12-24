''' moderator_commands.py A discord.py cog containing moderator commands for the Vyrtuous bot.

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
from collections import defaultdict
from typing import Optional
from discord import app_commands
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.discord_message_service import AppPaginator, DiscordMessageService, Paginator
from vyrtuous.utils.alias import Alias
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.duration import Duration
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.stage import Stage
from vyrtuous.utils.temporary_room import TemporaryRoom
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.voice_mute import VoiceMute
   
class ModeratorCommands(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.config = bot.config
        self.bot.db_pool = bot.db_pool
        self.emoji = Emojis()
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.channel_service = ChannelService()
        self.member_service = MemberService()
    
    # DONE
    @app_commands.command(name='bans', description='Lists ban statistics.')
    @app_commands.describe(target='"all", channel name/ID/mention, or user mention/ID')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_bans_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        member_obj = await self.member_service.resolve_member(interaction, target)
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content='\U0001F6AB Only owners or developers can list all bans across the server.')
            bans = await Ban.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not bans:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No active bans found in {interaction.guild.name}.')
            grouped = defaultdict(list)
            for ban in bans:
                grouped[ban.channel_snowflake].append(ban)
            pages = []
            for ch_id, records in grouped.items():
                ch = interaction.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Ban records for {ch_name}',
                    color=discord.Color.red()
                )
                for ban in bans:
                    user = interaction.guild.get_member(ban.member_snowflake)
                    reason = ban.reason or 'No reason provided'
                    if ban.expires_at is None:
                        duration_str = 'Permanent'
                    else:
                        now = discord.utils.utcnow()
                        delta = ban.expires_at - now
                        if delta.total_seconds() <= 0: duration_str = 'Expired'
                        else:
                            days, seconds = delta.days, delta.seconds
                            hours = seconds // 3600
                            minutes = (seconds % 3600) // 60
                            duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                    mention = user.name if user else f'`{record["discord_snowflake"]}`'
                    embed.add_field(name='User', value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        if member_obj:
            bans = await Ban.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            bans = [b for b in bans if interaction.guild.get_channel(b.channel_snowflake)]
            if not bans: return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not banned in any channels.')
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Ban records for {member_obj.display_name}',
                color=discord.Color.red()
            )
            for ban in bans:
                ch_obj = interaction.guild.get_channel(ban.channel_snowflake)
                channel_mention = ch_obj.mention if ch_obj else f'Channel ID `{ban.channel_snowflake}`'
                reason = ban.reason or 'No reason provided'
                if ban.expires_at is None:
                    duration_str = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = ban.expires_at - now
                    if delta.total_seconds() <= 0:
                        duration_str = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                embed.add_field(name=channel_mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
            return await interaction.response.send_message(embed=embed)
        elif channel_obj:
            bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not bans:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No active bans found for {channel_obj.mention}.')
            lines = []
            for ban in bans:
                uid = ban.member_snowflake
                member_obj = interaction.guild.get_member(uid)
                if not member_obj:
                    continue
                name = member_obj.display_name
                if ban.expires_at is None:
                    time_left = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = ban.expires_at - now
                    if delta.total_seconds() <= 0:
                        time_left = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        time_left = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                lines.append(f'• {name} — {time_left} — <@{uid}>')
            if not lines:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No active bans for users currently in {interaction.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Ban records for {channel_obj.name}',
                    description='\n'.join(chunk),
                    color=discord.Color.red()
                )
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        return await interaction.response.send_message(content='\U0001F6AB You must specify a member, a text channel or use "all".')
        
    # DONE
    @commands.command(name='bans', description='Lists ban statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_bans_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj and member_obj.id == ctx.guild.me.id:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot list bans on the bot.')
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        if member_obj:
            target = None
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list all bans across the server.')
            bans = await Ban.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not bans:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No active bans found in {ctx.guild.name}.')
            grouped = defaultdict(list)
            for ban in bans:
                grouped[ban.channel_snowflake].append(ban)
            pages = []
            for ch_id, records in grouped.items():
                ch = ctx.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                embed = discord.Embed(
                    title = f'{self.emoji.get_random_emoji()} Ban records for {ch_name}',
                    color = discord.Color.red()
                )
                for ban in bans:
                    user = ctx.guild.get_member(ban.member_snowflake)
                    reason = ban.reason or 'No reason provided'
                    if ban.expires_at is None:
                        duration_str='Permanent'
                    else:
                        now = discord.utils.utcnow()
                        delta = ban.expires_at - now
                        if delta.total_seconds() <= 0:
                            duration_str = 'Expired'
                        else:
                            days, seconds = delta.days,delta.seconds
                            hours = seconds // 3600
                            minutes = (seconds % 3600) // 60
                            duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                    mention = user.name if user else f'`{record["discord_snowflake"]}`'
                    embed.add_field(name='User', value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        elif member_obj:
            bans = await Ban.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if not bans:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not banned in any channels.')
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Ban records for {member_obj.display_name}',
                color=discord.Color.red()
            )
            for ban in bans:
                channel_obj = ctx.guild.get_channel(ban.channel_snowflake)
                channel_mention = channel_obj.mention if channel_obj else f'Channel ID `{ban.channel_snowflake}`'
                reason = ban.reason or 'No reason provided'
                if ban.expires_at is None:
                    duration_str = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = ban.expires_at - now
                    if delta.total_seconds() <= 0:
                        duration_str = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        duration_str = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                embed.add_field(name = channel_mention, value = f'Reason: {reason}\nDuration: {duration_str}', inline=False)
            return await self.handler.send_message(ctx, embed=embed)
        elif channel_obj:
            bans = await Ban.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not bans:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No active bans found for {channel_obj.mention}.')
            lines = []
            for ban in bans:
                uid = ban.member_snowflake
                member_obj = ctx.guild.get_member(uid)
                if not member_obj:
                    continue
                name = member_obj.display_name
                if ban.expires_at is None:
                    time_left = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = ban.expires_at - now
                    if delta.total_seconds() <= 0:
                        time_left = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        time_left = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                lines.append(f'• {name} — {time_left} — <@{uid}>')
            if not lines:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No active bans for users currently in {ctx.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i+chunk_size]
                embed = discord.Embed(
                    title = f'{self.emoji.get_random_emoji()} Ban records for {channel_obj.name}',
                    description = '\n'.join(chunk),
                    color = discord.Color.red()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify a member, a text channel or use "all".')
        
    
    # DONE
    @app_commands.command(name='caps', description='List active caps for a channel or all channels if "all" is provided.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    @app_commands.describe(target='"all", channel name/ID/mention')
    async def list_caps_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                return await interaction.response.send_message(content='\U0001F6AB Only owners or developers can list all caps.')
            caps = await Cap.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not caps:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No caps found server-wide.')
            lines = []
            for cap in caps:
                ch = interaction.guild.get_channel(cap.channel_snowflake)
                ch_name = ch.mention if ch else f'Channel ID `{cap.channel_snowflake}`'
                lines.append(f'**{cap.moderation_type} in {ch_name}** → `{Duration.convert_timedelta_seconds(cap.duration)}`')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                embed = discord.Embed(
                    title="{self.emoji.get_random_emoji()} All Active Caps in Server",
                    description="\n".join(lines[i:i+chunk_size]),
                    color=discord.Color.red()
                )
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
        if not caps:
            return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No caps found for this channel.')
        lines = []
        for cap in caps:
            lines.append(f'**{cap.moderation_type} in {channel_obj.mention}** → `{Duration.convert_timedelta_seconds(cap.duration)}`')
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} Active Caps for {channel_obj.mention}",
            description="\n".join(lines),
            color=discord.Color.red()
        )
        return await interaction.response.send_message(embed=embed)

    # DONE
    @commands.command(name='caps', help='List active caps for a channel or all channels if "all" is provided.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_caps_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention')
    ) -> None:
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list all caps.')
            caps = await Cap.fetch_by_guild(guild_snowflake=ctx.guild.id)
            if not caps:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No caps found server-wide.')
            lines = []
            for cap in caps:
                ch = ctx.guild.get_channel(cap.channel_snowflake)
                ch_name = ch.mention if ch else f'Channel ID `{cap.channel_snowflake}`'
                lines.append(f'**{cap.moderation_type} in {ch_name}** → `{Duration.convert_timedelta_seconds(cap.duration)}`')
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} All Active Caps in Server',
                description='\n'.join(lines),
                color=discord.Color.red()
            )
            return await self.handler.send_message(ctx, embed=embed)
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        caps = await Cap.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
        if not caps:
            return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No caps found for this channel.')
        lines = []
        for cap in caps:
            lines.append(f'**{cap.moderation_type} in {channel_obj.mention}** → `{Duration.convert_timedelta_seconds(cap.duration)}`')
        embed = discord.Embed(
            title=f'{self.emoji.get_random_emoji()} Active Caps for {channel_obj.mention}',
            description='\n'.join(lines),
            color=discord.Color.red()
        )
        await self.handler.send_message(ctx, embed=embed)
    
 
    # DONE
    @app_commands.command(name='cmds', description='List command aliases routed to a specific channel, temp room, or all channels if "all" is provided.')
    @app_commands.describe(target='"all", channel name/ID/mention')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_commands_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        found_aliases = False
        lines = []
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content='\U0001F6AB Only owners, developers or administrators can list all aliases across the server.')
            aliases = await Alias.fetch_by_guild(guild_snowflake=interaction.guild.id)
            if not aliases:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No aliases found.')
            lines.extend(Alias.format_aliases(aliases))
            found_aliases = len(lines) > 0
            if not found_aliases:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No aliases found in this server.')
            pages = []
            chunk_size = 18
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} All Aliases in Server', description='\n'.join(chunk), color=discord.Color.blue())
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        else:
            if channel_obj is None or channel_obj.type != discord.ChannelType.voice:
                return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=interaction.channel.id, guild_snowflake=interaction.guild.id)
            if not aliases:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No aliases found.')
            lines.extend(Alias.format_aliases(aliases))
            found_aliases = len(lines) > 0
            if not found_aliases:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No aliases found for the requested target: {target if target else channel_obj.mention}.')
        await interaction.response.send_message(embed=embed)
        
    # DONE
    @commands.command(name='cmds', help='List command aliases routed to a specific channel, temp room, or all channels if "all" is provided.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_commands_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or temp room name')
    ) -> None:
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        found_aliases = False
        lines = []
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners, developers and administrators can list all aliases across the server.')
            aliases = await Alias.fetch_by_guild(guild_snowflake=ctx.guild.id)
            lines.extend(Alias.format_aliases(aliases))
            found_aliases = len(lines) > 0
            if not found_aliases:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No aliases found in this server.')
            pages = []
            chunk_size = 18
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(title='All Aliases in Server', description='\n'.join(chunk), color=discord.Color.blue())
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        else:
            if channel_obj is None or channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            aliases = await Alias.fetch_by_channel_and_guild(channel_snowflake=ctx.channel.id, guild_snowflake=ctx.guild.id)
            if not aliases:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No aliases found.')
            lines.extend(Alias.format_aliases(aliases))
            found_aliases = len(lines) > 0
            if not found_aliases:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No aliases found for the requested target: {target if target else channel_obj.mention}.')
        embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} Aliases for {channel_obj.mention}', description='\n'.join(lines), color=discord.Color.blue())
        await self.handler.send_message(ctx, embed=embed)
 
    # DONE
    @app_commands.command(name='del', description='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @app_commands.describe(
        message='Message snowflake ID',
        channel='Tag a channel or include its snowflake ID'
    )
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def delete_message_app_command(
        self,
        interaction: discord.Interaction,
        message: AppMessageSnowflake,
        channel: AppChannelSnowflake
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        member_permission_role = await is_owner_developer_administrator_coordinator_moderator_via_channel_member(channel_obj, interaction.author)
        if member_permission_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator', 'Moderator'):
            return await interaction.response.send_message(content=f'\U0001F6AB You are not permitted to delete messages in {channel_obj.mention}.')
        try:
            msg = await channel_obj.fetch_message(message)
        except:
            logger.warning('No message with this ID exists.')
        if not msg:
            return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No message with ID `{message}` found in `{channel_obj.name}`.')
        try:
            await msg.delete()
            return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Message `{message}` deleted successfully.')
        except discord.Forbidden:
            logger.warning('Missing permissions to delete the message.')
        return await interaction.response.send_message(content='\U0001F6AB Failed to delete the message.')

    # DONE
    @commands.command(name='del', help='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message: MessageSnowflake = commands.parameter(default=None, description='Message snowflake'),
        *,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Channel or snowflake')
    ):
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        member_permission_role = await is_owner_developer_administrator_coordinator_moderator_via_channel_member(channel_obj, ctx.author)
        if member_permission_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator', 'Moderator'):
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not permitted to delete messages in {channel_obj.mention}.')
        try:
            msg = await channel_obj.fetch_message(message)
        except:
            logger.warning('No message with this ID exists.')
        if not msg:
            return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No message with ID `{message}` found in `{channel_obj.name}`.')
        try:
            await msg.delete()
            return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Message `{message}` deleted successfully.')
        except discord.Forbidden:
            logger.warning('Missing permissions to delete the message.')
        return await self.handler.send_message(ctx, content='\U0001F6AB Failed to delete the message.')
        
    # DONE
    @app_commands.command(name='flags', description='List flag statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_flags_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        member_obj = await self.member_service.resolve_member(interaction, target)
        if member_obj and member_obj.id == interaction.guild.me.id:
            return await interaction.response.send_message(content='\U0001F6AB You cannot list flags on the bot.')
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        if member_obj:
            target=None
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content='\U0001F6AB Only owners, developers or administrators can list flags across all channels.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('SELECT channel_id, discord_snowflake FROM active_flags WHERE guild_id=$1', interaction.guild.id)
            if not rows: return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No users are flagged in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows: channel_map[row['channel_id']].append(row['discord_snowflake'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                ch = interaction.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Flag records for {ch_name}',
                    color=discord.Color.yellow()
                )
                for uid in user_ids:
                    m = interaction.guild.get_member(uid)
                    mention = m.mention if m else f'<@{uid}>'
                    embed.add_field(name=f'{interaction.guild.name}', value=f'• {mention}', inline=False)
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('SELECT channel_id, reason FROM active_flags WHERE guild_id=$1 AND discord_snowflake=$2', interaction.guild.id, member_obj.id)
            rows = [r for r in rows if interaction.guild.get_channel(r['channel_id'])]
            if not rows:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not flagged in any voice channels.')
            lines = []
            for r in rows:
                ch = interaction.guild.get_channel(r['channel_id'])
                ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                reason = r['reason'] or "No reason given"
                lines.append(f'• {ch_name} — {reason}')
            embed=discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Flag records for {member_obj.display_name}',
                description='\n'.join(lines),
                color=discord.Color.orange()
            )
            return await interaction.response.send_message(embed=embed)
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('SELECT discord_snowflake FROM active_flags WHERE guild_id=$1 AND channel_id=$2', interaction.guild.id, channel_obj.id)
            if not rows:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No users are flagged for {channel_obj.mention}.')
            pages = []
            chunk_size = 18
            for i in range(0, len(rows), chunk_size):
                chunk = rows[i:i+chunk_size]
                formatted_lines = []
                for record in chunk:
                    uid = record['discord_snowflake']
                    m = interaction.guild.get_member(uid)
                    formatted_lines.append(f'• {m.display_name} — <@{uid}>')
                if formatted_lines:
                    embed=discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Flag records for {channel_obj.name}',
                        color=discord.Color.red()
                    )
                    embed.add_field(name=f'{interaction.guild.name}', value='\n'.join(formatted_lines), inline=False)
                    pages.append(embed)
            if not pages:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No flagged users currently in {interaction.guild.name}.')
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
    
    # DONE
    @commands.command(name='flags', help='List flag statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_flags_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj and member_obj.id == ctx.guild.me.id:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot list flags on the bot.')
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        if member_obj:
            target = None
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners, developers or administrators can list flags across all channels.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT channel_id, discord_snowflake
                    FROM active_flags
                    WHERE guild_id = $1
                ''', ctx.guild.id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No users are flagged in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows:
                channel_map[row['channel_id']].append(row['discord_snowflake'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                ch = ctx.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Flag records for {ch_name}',
                    color=discord.Color.yellow()
                )
                for uid in user_ids:
                    m = ctx.guild.get_member(uid)
                    mention = m.mention if m else f'<@{uid}>'
                    embed.add_field(name=f'{ctx.guild.name}', value=f'• {mention}', inline=False)
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT channel_id, reason
                    FROM active_flags
                    WHERE guild_id = $1 AND discord_snowflake = $2
                ''', ctx.guild.id, member_obj.id)
                rows = [r for r in rows if ctx.guild.get_channel(r['channel_id'])]
                if not rows:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not flagged in any voice channels.')
                lines = []
                for r in rows:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    reason = r['reason'] or "No reason given"
                    lines.append(f'• {ch_name} — {reason}')
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Flag records for {member_obj.display_name}',
                    description='\n'.join(lines),
                    color=discord.Color.orange()
                )
                return await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.all())
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT discord_snowflake
                    FROM active_flags
                    WHERE guild_id = $1 AND channel_id = $2
                ''', ctx.guild.id, channel_obj.id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No users are flagged for {channel_obj.mention}.')
            pages = []
            chunk_size = 18
            for i in range(0, len(rows), chunk_size):
                chunk = rows[i:i+chunk_size]
                formatted_lines = []
                for record in chunk:
                    uid = record['discord_snowflake']
                    member_obj = ctx.guild.get_member(uid)
                    formatted_lines.append(f'• {member_obj.display_name} — <@{uid}>')
                if formatted_lines:
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Flag records for {channel_obj.name}',
                        color=discord.Color.red()
                    )
                    embed.add_field(name=f'{ctx.guild.name}', value='\n'.join(formatted_lines), inline=False)
                    pages.append(embed)
            if not pages:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No flagged users currently in {ctx.guild.name}.')
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        
    # DONE
    @app_commands.command(name='ls', description='List users cowed as going vegan in this guild.')
    @app_commands.describe(target='A member or channel snowflake ID/mention')
    async def list_members_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        member_obj = await self.member_service.resolve_member(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        async with self.bot.db_pool.acquire() as conn:
            if member_obj:
                rows = await conn.fetch('''SELECT channel_id, created_at FROM active_cows WHERE guild_id=$1 AND discord_snowflake=$2''', interaction.guild.id, member_obj.id)
                if not rows:
                    return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not cowed in any channels.')
                lines = []
                for r in rows:
                    ch = interaction.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    created_at = discord.utils.format_dt(r['created_at'],style='R') if r['created_at'] else ''
                    lines.append(f'• {ch_name} — {created_at}')
                embed=discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} {member_obj.display_name}',
                    description='\n'.join(lines),
                    color=discord.Color.green()
                )
                return await interaction.response.send_message(embed=embed, allowed_mentions=discord.AllowedMentions.all())
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
                rows = await conn.fetch('''SELECT discord_snowflake, created_at FROM active_cows WHERE guild_id=$1 AND channel_id=$2''', interaction.guild.id, channel_obj.id)
                if not rows:
                    return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No users are cowed in {channel_obj.mention}.')
                lines = []
                for row in rows:
                    uid = row['discord_snowflake']
                    m = interaction.guild.get_member(uid)
                    if not m:
                        continue
                    created_at = discord.utils.format_dt(row['created_at'],style='R') if row['created_at'] else ''
                    lines.append(f'• {m.display_name} — <@{uid}> — {created_at}')
                if not lines:
                    return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No new vegans currently in {interaction.guild.name}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk=lines[i:i+chunk_size]
                    embed=discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Vegan records for {channel_obj.name}',
                        description='\n'.join(chunk),
                        color=discord.Color.green()
                    )
                    pages.append(embed)
                paginator = AppPaginator(self.bot, interaction, pages)
                return await paginator.start()
                
    # DONE
    @commands.command(name='ls', help='List users cowed as going vegan in this guild.')
    async def list_members_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj:
            target = None
        if member_obj and member_obj.id == ctx.guild.me.id:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot get a vegan status for the bot.')
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        async with self.bot.db_pool.acquire() as conn:
            if member_obj:
                rows = await conn.fetch('''
                    SELECT channel_id, created_at
                    FROM active_cows
                    WHERE guild_id = $1 AND discord_snowflake = $2
                ''', ctx.guild.id, member_obj.id)
                if not rows:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not cowed in any channels.')
                lines = []
                for r in rows:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    created_at = discord.utils.format_dt(r['created_at'], style='R') if r['created_at'] else ''
                    lines.append(f'• {ch_name} — {created_at}')
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} {member_obj.display_name}',
                    description='\n'.join(lines),
                    color=discord.Color.green()
                )
                return await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.all())
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
                rows = await conn.fetch('''
                    SELECT discord_snowflake, created_at
                    FROM active_cows
                    WHERE guild_id = $1 AND channel_id = $2
                ''', ctx.guild.id, channel_obj.id)
                if not rows:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No users are cowed in {channel_obj.mention}.')
                lines = []
                for row in rows:
                    uid = row['discord_snowflake']
                    m = ctx.guild.get_member(uid)
                    if not m:
                        continue
                    created_at = discord.utils.format_dt(row['created_at'], style='R') if row['created_at'] else ''
                    lines.append(f'• {m.display_name} — <@{uid}> — {created_at}')
                if not lines:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No new vegans currently in {ctx.guild.name}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i + chunk_size]
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Vegan records for {channel_obj.name}',
                        description='\n'.join(chunk),
                        color=discord.Color.green()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()

 # DONE
    @app_commands.command(name='migrate', description='Migrate a temporary room to a new channel.')
    @app_commands.describe(
        old_name='Old temporary room name',
        channel='New channel to migrate to'
    )
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def migrate_temp_room_app_command(
        self,
        interaction: discord.Interaction,
        old_name: str,
        channel: AppChannelSnowflake
    ):
        old_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=interaction.guild.id, room_name=old_name)
        if old_room:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
            is_owner = old_room.member_snowflake == interaction.user.id
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
            if highest_role not in ('Owner', 'Developer', 'Administrator') or not is_owner:
                return await interaction.response.send_message(content=f'You are not the owner nor have elevated permission to migrate rooms.')
            await TemporaryRoom.update_by_source_and_target(guild_snowflake=interaction.guild.id, room_name=channel_obj.id, source_channel_snowflake=old_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            new_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=interaction.guild.id, room_name=channel_obj.name)
            await Alias.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Ban.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Cap.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Coordinator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Moderator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Stage.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await TextMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await VoiceMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
        return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Temporary room {old_name} migrated to {channel}.')
    
    # DONE
    @commands.command(name='migrate', help='Migrate a temporary room to a new channel by snowflake.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def migrate_temp_room_text_command(
        self,
        ctx: commands.Context,
        old_name: Optional[str] = commands.parameter(default=None, description='Provide a channel name'),
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        old_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=ctx.guild.id, room_name=old_name)
        if old_room:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            if highest_role not in ('Owner', 'Developer', 'Administrator') or not is_owner:
                return await self.handler.send_message(ctx, content=f'You are not the owner nor have elevated permission to migrate rooms.')
            await TemporaryRoom.update_by_source_and_target(guild_snowflake=ctx.guild.id, room_name=channel_obj.name, source_channel_snowflake=old_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            new_room = await TemporaryRoom.fetch_by_guild_and_room_name(guild_snowflake=ctx.guild.id, room_name=channel_obj.name)
            await Alias.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Ban.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Cap.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Coordinator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Moderator.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await Stage.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await TextMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
            await VoiceMute.update_by_source_and_target(source_channel_snowflake=new_room.channel_snowflake, target_channel_snowflake=channel_obj.id)
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Temporary room `{old_name}` migrated to {channel}.')

    # DONE
    @app_commands.command(name='mutes', description='Lists mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        member_obj = await self.member_service.resolve_member(interaction, target)
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content=f'\U0001F6AB You do not have permission to use this command (`mutes`) in {channel_obj.mention}.')
            voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=interaction.guild.id, target="user")
            if not voice_mutes:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No muted users currently in {interaction.guild.name}.')
            grouped = defaultdict(list)
            for voice_mute in voice_mutes:
                grouped[voice_mute.channel_snowflake].append(voice_mute)
            pages = []
            for channel_id, user_entries in sorted(grouped.items()):
                channel = interaction.guild.get_channel(channel_id)
                channel_name = channel.mention if channel else f'Unknown Channel ({channel_id})'
                chunk_size = 18
                for i in range(0, len(user_entries), chunk_size):
                    embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} Mutes records for {channel_name}', color=discord.Color.orange())
                    for voice_mute in user_entries[i:i + chunk_size]:
                        user_id = voice_mute.member_snowflake
                        member = interaction.guild.get_member(user_id)
                        name = member.display_name if member else f'User ID {user_id}'
                        mention = member.mention if member else f'`{user_id}`'
                        duration_str = Duration.output_display_from_datetime(r['expires_at'])
                        embed.add_field(name=name, value=f'{mention}\nReason: {r["reason"]}\nDuration: {duration_str}', inline=False)
                    pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        elif member_obj:
            voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, target="user")
            voice_mutes = [voice_mute for voice_mute in voice_mutes if interaction.guild.get_channel(voice_mute.channel_snowflake)]
            if not voice_mutes:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not muted in any voice channels.')
            lines = []
            for voice_mute in voice_mutes:
                ch = interaction.guild.get_channel(voice_mute.channel_snowflake)
                ch_mention = ch.mention if ch else f'`{voice_mute.channel_snowflake}`'
                lines.append(f'• {ch_mention} — {r["reason"]} — {Duration.output_display_from_datetime(voice_mute.expires_at)}')
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Mute records for {member_obj.display_name}',
                description='\n'.join(lines),
                color=discord.Color.orange()
            )
            return await interaction.response.send_message(embed=embed)
        elif channel_obj:
            voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target="user")
            if not voice_mutes:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No users are currently muted in {channel_obj.mention}.')
            lines = []
            for voice_mute in voice_mutes:
                uid = voice_mute.member_snowflake
                member_obj = interaction.guild.get_member(uid)
                if not member_obj:
                    continue
                lines.append(f'• {member_obj.display_name} — <@{uid}> — {Duration.output_display_from_datetime(voice_mute.expires_at)}')
            if not lines:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No muted users currently in {interaction.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Mute records for {channel_obj.name}',
                    color=discord.Color.orange()
                )
                embed.add_field(name=f'{interaction.guild.name}', value='\n'.join(chunk), inline=False)
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()

    # DONE
    @commands.command(name='mutes', help='Lists mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`mutes`) in {channel_obj.mention}.')
            voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=ctx.guild.id, target="user")
            if not voice_mutes:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No muted users currently in {ctx.guild.name}.')
            grouped = defaultdict(list)
            for voice_mute in voice_mutes:
                grouped[voice_mute.channel_snowflake].append(voice_mute)
            pages = []
            for channel_id, user_entries in sorted(grouped.items()):
                channel = ctx.guild.get_channel(channel_id)
                channel_name = channel.mention if channel else f'Unknown Channel ({channel_id})'
                chunk_size = 18
                for i in range(0, len(user_entries), chunk_size):
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Mutes records for {channel_name}',
                        color=discord.Color.orange()
                    )
                    for voice_mute in user_entries[i:i + chunk_size]:
                        user_id = voice_mute.member_snowflake
                        member = ctx.guild.get_member(user_id)
                        name = member.display_name if member else f'User ID {user_id}'
                        mention = member.mention if member else f'`{user_id}`'
                        reason = record['reason'] or 'No reason provided'
                        duration_str = Duration.output_display_from_datetime(voice_mute.expires_at)
                        embed.add_field(name=name, value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                    pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        if member_obj:
            voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, target="user")
            if not voice_mutes:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not muted in any voice channels.')
            description_lines = []
            for voice_mute in voice_mutes:
                channel_obj = ctx.guild.get_channel(voice_mute.channel_snowflake)
                channel_mention = channel_obj.mention if channel_obj else f'`{voice_mute.channel_snowflake}`'
                reason = record['reason']
                duration_str = Duration.output_display_from_datetime(voice_mute.expires_at)
                description_lines.append(f'• {channel_mention} — {reason} — {duration_str}')
            embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} Mute records for {member_obj.display_name}', description='\n'.join(description_lines), color=discord.Color.orange())
            return await self.handler.send_message(ctx, embed=embed)
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target="user")
                if not voice_mutes:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No users are currently muted in {channel_obj.mention}.')
                description_lines = []
                for voice_mute in voice_mutes:
                    uid = voice_mute.member_snowflake
                    member_obj = ctx.guild.get_member(uid)
                    if not member_obj:
                        continue
                    name = member_obj.display_name
                    duration_str = Duration.output_display_from_datetime(voice_mute.expires_at)
                    description_lines.append(f'• {name} — <@{uid}> — {duration_str}')
                if not description_lines:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No muted users currently in {channel_obj.mention}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(description_lines), chunk_size):
                    chunk = description_lines[i:i + chunk_size]
                    embed = discord.Embed(title=f'{self.emoji.get_random_emoji()} Mute records for {channel_obj.name}', color=discord.Color.orange())
                    embed.add_field(name=f'{ctx.guild.name}', value='\n'.join(chunk), inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
    
    # DONE
    @app_commands.command(name='mstage', description='Mute/unmute a member in the active stage.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    async def stage_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake
    ):
        member_obj = await self.member_service.resolve_member(interaction, member)
        if member_obj:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not stage:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No active stage found.')
            success = await has_equal_or_higher_role(interaction, member_obj, channel_obj)
            if not success:
                return await interaction.response.send_message(content=f"\U0001F6AB You are not allowed to mute/unmute {member_obj.mention} because they are a higher/or equivalent role than you in {channel_obj.mention}.")
            try:
                await member_obj.edit(mute=not member_obj.voice.mute)
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been {"muted" if member_obj.voice.mute else "unmuted"}.')
            except Exception as e:
                logger.warning(f'Failed to toggle mute: {e}')
        return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Toggled mute.')
                
    # DONE
    @commands.command(name='mstage', help='Mute/unmute a member in the active stage.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def stage_mute_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: ChannelSnowflake = commands.parameter(default=None, description="Tag a channel or include it's snowflake ID")
    ) -> None:
        member_obj = await self.member_service.resolve_member(ctx, member)
        if member_obj:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not stage:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No active stage found.')
            success = await has_equal_or_higher_role(ctx, member_obj, channel_obj)
            if not success:
                return await self.handler.send_message(ctx, content=f"\U0001F6AB You are not allowed to mute/unmute {member_obj.mention} because they are a higher/or equivalent role than you in {channel_obj.mention}.")
            try:
                await member_obj.edit(mute=not member_obj.voice.mute)
            except Exception as e:
                logger.warning(f'Failed to toggle mute: {e}')
        return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Toggled mute.')
    
    # DONE
    @app_commands.command(name='stages', description='Lists stage mute statistics.')
    @app_commands.describe(target='"all", channel name/ID/mention')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_stages_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await interaction.response.send_message(content='\U0001F6AB Only owners/devs can view all stages.')
                stages = await Stage.fetch_by_guild(guild_snowflake=interaction.guild.id)
                if not stages:
                    return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No active stages in {interaction.guild.name}.')
                pages, chunk_size = [], 8
                for i in range(0, len(stages), chunk_size):
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Active Stages in {interaction.guild.name}',
                        color=discord.Color.purple()
                    )
                    for s in stages[i:i+chunk_size]:
                        ch = interaction.guild.get_channel(s.channel_snowflake)
                        ch_name = ch.mention if ch else f'Unknown Channel ({s.channel_snowflake})'
                        voice_mutes = await VoiceMute.fetch_by_guild_and_target(guild_snowflake=interaction.guild.id, target="room")
                        for voice_mute in voice_mutes:
                            embed.add_field(name=ch_name, value=f'Active stage mutes: {voice_mute.member_snowflake}', inline=False)
                    pages.append(embed)
                paginator = AppPaginator(self.bot, interaction, pages)
                return await paginator.start()
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
            if not stage:
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No active stage in {channel_obj.mention}.')
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target="room")
            initiator = interaction.guild.get_member(stage.member_snowflake)
            initiator_name = initiator.mention if initiator else f'`{stage.member_snowflake}`'
            expires = Duration.output_display_from_datetime(stage.expires_at) if stage.expires_at else 'No expiration'
            lines = []
            for m in voice_mutes:
                user = interaction.guild.get_member(m.member_snowflake)
                duration_str = Duration.output_display_from_datetime(m.expires_at) if m.expires_at else 'No expiration'
                reason = m.reason or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = f'Initiator: {initiator_name}\nStage Expires: {expires}\nActive stage mutes: {len(lines)}\n\n'+'\n'.join(lines)
            pages, chunk_size = [], 18
            for i in range(0, len(description.splitlines()), chunk_size):
                embed=discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Stage info for {channel_obj.mention}',
                    description='\n'.join(description.splitlines()[i:i+chunk_size]),
                    color=discord.Color.purple()
                )
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
            
    # DONE
    @commands.command(name='stages', help='Lists stage mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_stages_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower() == 'all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Only owners/devs can view all stages.')
                stages = await Stage.fetch_by_guild(guild_snowflake=ctx.guild.id)
                if not stages:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No active stages in {ctx.guild.name}.')
                pages, chunk_size = [], 8
                for i in range(0, len(stages), chunk_size):
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Active Stages in {ctx.guild.name}',
                        color=discord.Color.purple()
                    )
                    voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target="room")
                    for s in stages[i:i+chunk_size]:
                        ch = ctx.guild.get_channel(s.channel_snowflake)
                        ch_name = ch.mention if ch else f'Unknown Channel ({s.channel_snowflake})'
                        for voice_mute in voice_mutes:
                            embed.add_field(name=ch_name, value=f'Active stage mutes: {voice_mute.member_snowflake}', inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            stage = await Stage.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
            if not stage:
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No active stage in {channel_obj.mention}.')
            voice_mutes = await VoiceMute.fetch_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target="room")
            initiator = ctx.guild.get_member(stage.member_snowflake)
            initiator_name = initiator.mention if initiator else f'`{stage.member_snowflake}`'
            expires = Duration.output_display_from_datetime(stage.expires_at) if stage.expires_at else 'No expiration'
            lines = []
            for voice_mute in voice_mutes:
                user = ctx.guild.get_member(voice_mute.member_snowflake)
                duration_str = Duration.output_display_from_datetime(voice_mute.expires_at) if voice_mute.expires_at else 'No expiration'
                reason = voice_mute.reason or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = (
                f'Initiator: {initiator_name}\n'
                f'Stage Expires: {expires}\n'
                f'Active stage mutes: {len(lines)}\n\n' + '\n'.join(lines)
            )
            pages, chunk_size = [], 18
            for i in range(0, len(description.splitlines()), chunk_size):
                embed = discord.Embed(
                    title=f'{self.emoji.get_random_emoji()} Stage info for {channel_obj.mention}',
                    description='\n'.join(description.splitlines()[i:i+chunk_size]),
                    color=discord.Color.purple()
                )
                pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
            
    # DONE
    @app_commands.command(name='tmutes', description='Lists text-mute statistics.')
    @app_commands.describe(target='"all", channel name/ID/mention, or user mention/ID')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_text_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        member_obj = await self.member_service.resolve_member(interaction, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower()=='all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await interaction.response.send_message(content='\U0001F6AB Only owners, developers and administrators can list all text-mutes across the server.')
                text_mutes = await TextMute.fetch_by_guild(guild_snowflake=interaction.guild.id)
                if not text_mutes:
                    return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No users are currently text-muted in {interaction.guild.name}.')
                grouped = defaultdict(list)
                for text_mute in text_mutes:
                    grouped[text_mute.channel_snowflake].append(text_mute)
                pages, chunk_size = [], 18
                for ch_id, entries in sorted(grouped.items()):
                    ch = interaction.guild.get_channel(ch_id)
                    ch_name = ch.mention
                    for i in range(0, len(entries), chunk_size):
                        embed= discord.Embed(
                            title=f'{self.emoji.get_random_emoji()} Text-mute records for {ch_name}',
                            color=discord.Color.orange()
                        )
                        for e in entries[i:i+chunk_size]:
                            user = interaction.guild.get_member(e.member_snowflake)
                            mention = user.name if user else f'`{e.member_snowflake}`'
                            reason = e.reason or 'No reason provided'
                            duration_str = Duration.output_display_from_datetime(e.expires_at)
                            embed.add_field(name=mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
                        pages.append(embed)
                paginator = AppPaginator(self.bot, interaction, pages)
                return await paginator.start()
            elif member_obj:
                text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                if not text_mutes:
                    return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} is not text-muted in any channels.')
                lines = []
                for text_mute in text_mutes:
                    ch = interaction.guild.get_channel(text_mutes.channel_snowflake)
                    duration_str = Duration.output_display_from_datetime(text_mute.expires_at)
                    lines.append(f'• {ch.mention} — {text_mutes.reason} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Text-mute records for {member_obj.display_name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = AppPaginator(self.bot, interaction, pages)
                return await paginator.start()
            elif channel_obj:
                text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id)
                if not text_mutes:
                    return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} No users are currently text-muted in {channel_obj.mention}.')
                lines = []
                for text_mute in text_mutes:
                    user = interaction.guild.get_member(text_mute.member_snowflake)
                    if not user:
                        continue
                    duration_str = Duration.output_display_from_datetime(text_mute.expires_at)
                    lines.append(f'• {user.mention} — {text_mutes.reason} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Text-mute records for {channel_obj.name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = AppPaginator(self.bot, interaction, pages)
                return await paginator.start()
        return await interaction.response.send_message(content='\U0001F6AB You must specify "all", a member, or a text channel.')
    
    # DONE
    @commands.command(name='tmutes', help='Lists text-mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_text_mutes_text_command(
        self,
        ctx: commands.Context,
        target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or user mention/ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        async with self.bot.db_pool.acquire() as conn:
            if target and target.lower()=='all':
                if highest_role not in ('Owner', 'Developer', 'Administrator'):
                    return await self.handler.send_message(ctx, content='\U0001F6AB Only owners, developers and administrators can list all text-mutes across the server.')
                text_mutes = await TextMute.fetch_by_guild(guild_snowflake=ctx.guild.id)
                if not text_mutes:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No users are currently text-muted in {ctx.guild.name}.')
                grouped = defaultdict(list)
                for text_mute in text_mutes:
                    grouped[text_mute.channel_snowflake].append(text_mute)
                pages, chunk_size = [], 18
                for ch_id, entries in sorted(grouped.items()):
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    for i in range(0, len(entries), chunk_size):
                        embed = discord.Embed(
                            title=f'{self.emoji.get_random_emoji()} Text-mute records for {ch_name}',
                            color=discord.Color.orange()
                        )
                        for e in entries[i:i + chunk_size]:
                            user = ctx.guild.get_member(e.discord_snowflake)
                            mention = user.name if user else f'<@{e.discord_snowflake}>'
                            reason = getattr(e, 'reason', 'No reason provided')
                            expires_at = getattr(e, 'expires_at', None)
                            duration_str = Duration.output_display_from_datetime(expires_at) if expires_at else "Permanent"
                            embed.add_field(
                                name=mention,
                                value=f'Reason: {reason}\nDuration: {duration_str}',
                                inline=False
                            )
                        pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif member_obj:
                text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                if not text_mutes:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not text-muted in any channels.')
                lines = []
                for r in text_mutes:
                    ch = ctx.guild.get_channel(r.channel_snowflake)
                    ch_name = ch.mention if ch else f'Channel ID `{r.channel_snowflake}`'
                    duration_str = Duration.output_display_from_datetime(r.expires_at) if r.expires_at else "Permanent"
                    reason = r.reason if r.reason else 'No reason provided'
                    lines.append(f'• {ch_name} — {reason} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Text-mute records for {member_obj.display_name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif channel_obj:
                text_mutes = await TextMute.fetch_by_channel_and_guild(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id)
                if not text_mutes:
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} No users are currently text-muted in {channel_obj.mention}.')
                lines = []
                for r in text_mutes:
                    user = ctx.guild.get_member(r.discord_snowflake)
                    mention = user.mention if user else f'<@{r.discord_snowflake}>'
                    duration_str = Duration.output_display_from_datetime(r.expires_at) if r.expires_at else "Permanent"
                    reason = r.reason if r.reason else 'No reason provided'
                    lines.append(f'• {mention} — {reason} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'{self.emoji.get_random_emoji()} Text-mute records for {channel_obj.name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
        return await self.handler.send_message(ctx, content='\U0001F6AB You must specify "all", a member, or a text channel.')
        
async def setup(bot: DiscordBot):
    cog = ModeratorCommands(bot)
    await bot.add_cog(cog)

