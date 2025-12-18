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
from vyrtuous.utils.cap import Cap
from vyrtuous.utils.duration import Duration
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.setup_logging import logger
   
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
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(interaction, target)
        if member_obj and member_obj.id == interaction.guild.me.id:
            return await interaction.response.send_message(content='\U0001F6AB You cannot list bans on the bot.')
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        if member_obj:
            target = None
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content='\U0001F6AB Only owners or developers can list all bans across the server.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''
                    SELECT discord_snowflake, room_name, channel_id, expires_at, reason
                    FROM active_bans
                    WHERE guild_id=$1
                    ORDER BY channel_id, room_name, expires_at NULLS LAST
                ''', interaction.guild.id)
            if not rows:
                return await interaction.response.send_message(content=f'\U0001F6AB No active bans found in {interaction.guild.name}.')
            grouped = defaultdict(list)
            for row in rows: grouped[row['channel_id']].append(row)
            pages = []
            for ch_id, records in grouped.items():
                ch = interaction.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                embed = discord.Embed(
                    title=f'\u26D4 Ban records for {ch_name}',
                    color=discord.Color.red()
                )
                for record in records:
                    user = interaction.guild.get_member(record['discord_snowflake'])
                    reason = record['reason'] or 'No reason provided'
                    if record['expires_at'] is None:
                        duration_str = 'Permanent'
                    else:
                        now = discord.utils.utcnow()
                        delta = record['expires_at'] - now
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
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''
                    SELECT channel_id, expires_at, reason
                    FROM active_bans
                    WHERE guild_id=$1 AND discord_snowflake=$2 AND room_name = $3
                ''', interaction.guild.id, member_obj.id, channel_obj.name)
            bans = [b for b in bans if interaction.guild.get_channel(b['channel_id'])]
            if not bans: return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} is not banned in any channels.', allowed_mentions=discord.AllowedMentions.none())
            embed = discord.Embed(
                title=f'\u26D4 Ban records for {member_obj.display_name}',
                color=discord.Color.red()
            )
            for record in bans:
                ch_obj = interaction.guild.get_channel(record['channel_id'])
                channel_mention = ch_obj.mention if ch_obj else f'Channel ID `{record["channel_id"]}`'
                reason = record['reason'] or 'No reason provided'
                if record['expires_at'] is None:
                    duration_str = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = record['expires_at'] - now
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
            if channel_obj.type != discord.ChannelType.voice:
                return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''
                    SELECT discord_snowflake, expires_at, reason
                    FROM active_bans
                    WHERE guild_id=$1 AND channel_id=$2 AND room_name = $3
                    ORDER BY expires_at NULLS LAST
                ''', interaction.guild.id, channel_obj.id, channel_obj.name)
            if not bans:
                return await interaction.response.send_message(content=f'\U0001F6AB No active bans found for {channel_obj.mention}.')
            lines = []
            for record in bans:
                uid = record['discord_snowflake']
                member_obj = interaction.guild.get_member(uid)
                if not member_obj: continue
                name = member_obj.display_name
                if record['expires_at'] is None:
                    time_left = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = record['expires_at'] - now
                    if delta.total_seconds() <= 0:
                        time_left = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        time_left = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                lines.append(f'• {name} — {time_left} — <@{uid}>')
            if not lines:
                return await interaction.response.send_message(content=f'\U0001F6AB No active bans for users currently in {interaction.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'\u26D4 Ban records for {channel_obj.name}',
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
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
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
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('''SELECT discord_snowflake, channel_id, room_name, expires_at, reason FROM active_bans WHERE guild_id = $1 ORDER BY channel_id, room_name, expires_at NULLS LAST''', ctx.guild.id)
            if not rows:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No active bans found in {ctx.guild.name}.')
            grouped = defaultdict(list)
            for row in rows:
                grouped[row['channel_id']].append(row)
            pages = []
            for ch_id, records in grouped.items():
                ch = ctx.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                embed = discord.Embed(
                    title = f'\u26D4 Ban records for {ch_name}',
                    color = discord.Color.red()
                )
                for record in records:
                    user = ctx.guild.get_member(record['discord_snowflake'])
                    reason = record['reason'] or 'No reason provided'
                    if record['expires_at'] is None:
                        duration_str='Permanent'
                    else:
                        now = discord.utils.utcnow()
                        delta = record['expires_at'] - now
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
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''SELECT channel_id,expires_at,reason FROM active_bans WHERE guild_id=$1 AND discord_snowflake=$2 AND room_name = $3''', ctx.guild.id, member_obj.id, channel_obj.name)
            bans = [b for b in bans if ctx.guild.get_channel(b['channel_id'])]
            if not bans:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not banned in any channels.', allowed_mentions=discord.AllowedMentions.none())
            embed = discord.Embed(
                title=f'\u26D4 Ban records for {member_obj.display_name}',
                color=discord.Color.red()
            )
            for record in bans:
                channel_obj = ctx.guild.get_channel(record['channel_id'])
                channel_mention = channel_obj.mention if channel_obj else f'Channel ID `{record["channel_id"]}`'
                reason = record['reason'] or 'No reason provided'
                if record['expires_at'] is None:
                    duration_str = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = record['expires_at'] - now
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
            if channel_obj.type != discord.ChannelType.voice:
               return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                bans = await conn.fetch('''SELECT discord_snowflake,expires_at,reason FROM active_bans WHERE guild_id=$1 AND channel_id=$2 AND room_name = $3 ORDER BY expires_at NULLS LAST''', ctx.guild.id, channel_obj.id, channel_obj.name)
            if not bans:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No active bans found for {channel_obj.mention}.')
            lines = []
            for record in bans:
                uid = record['discord_snowflake']
                member_obj = ctx.guild.get_member(uid)
                if not member_obj:
                    continue
                name = member_obj.display_name
                if record['expires_at'] is None:
                    time_left = 'Permanent'
                else:
                    now = discord.utils.utcnow()
                    delta = record['expires_at'] - now
                    if delta.total_seconds() <= 0:
                        time_left = 'Expired'
                    else:
                        days, seconds = delta.days, delta.seconds
                        hours = seconds // 3600
                        minutes = (seconds % 3600) // 60
                        time_left = f'{days}d {hours}h left' if days > 0 else f'{hours}h {minutes}m left' if hours > 0 else f'{minutes}m left'
                lines.append(f'• {name} — {time_left} — <@{uid}>')
            if not lines:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No active bans for users currently in {ctx.guild.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i+chunk_size]
                embed = discord.Embed(
                    title = f'\u26D4 Ban records for {channel_obj.name}',
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
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer'):
                return await interaction.response.send_message(content='\U0001F6AB Only owners or developers can list all caps.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(
                    "SELECT channel_id, duration_seconds, moderation_type FROM active_caps WHERE guild_id=$1",
                    interaction.guild.id
                )
            if not rows:
                return await interaction.response.send_message(content='\U0001F6AB No caps found server-wide.')
            lines = []
            for row in rows:
                ch = interaction.guild.get_channel(row["channel_id"])
                ch_name = ch.mention if ch else f'Channel ID `{row["channel_id"]}`'
                lines.append(f'**{row["moderation_type"]} in {ch_name}** → `{Duration.convert_timedelta_seconds(row["duration_seconds"])}`')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                embed = discord.Embed(
                    title="All Active Caps in Server",
                    description="\n".join(lines[i:i+chunk_size]),
                    color=discord.Color.red()
                )
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        caps = await Cap.get_caps_for_channel(interaction.guild.id, channel_obj.id)
        if not caps:
            return await interaction.response.send_message(content='\U0001F6AB No caps found for this channel.')
        lines = [
            f'**{moderation_type} in {channel_obj.mention}** → `{Duration.convert_timedelta_seconds(duration_seconds)}`'
            for duration_seconds, moderation_type in caps
        ]
        embed = discord.Embed(
            title=f"Active Caps for {channel_obj.mention}",
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
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners or developers can list all caps.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch(
                    "SELECT channel_id, duration_seconds, moderation_type FROM active_caps WHERE guild_id=$1",
                    ctx.guild.id
                )
            if not rows:
                return await self.handler.send_message(ctx, content='\U0001F6AB No caps found server-wide.')
            lines = []
            for row in rows:
                ch = ctx.guild.get_channel(row["channel_id"])
                ch_name = ch.mention if ch else f'Channel ID `{row["channel_id"]}`'
                lines.append(f'**{row["moderation_type"]} in {ch_name}** → `{Duration.convert_timedelta_seconds(row["duration_seconds"])}`')
            embed = discord.Embed(
                title='All Active Caps in Server',
                description='\n'.join(lines),
                color=discord.Color.red()
            )
            return await self.handler.send_message(ctx, embed=embed)
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        caps = await Cap.get_caps_for_channel(ctx.guild.id, channel_obj.id)
        if not caps:
            return await self.handler.send_message(ctx, content='\U0001F6AB No caps found for this channel.')
        lines = [f'**{mtype} in {channel_obj.mention}** → `{Duration.convert_timedelta_seconds(duration_seconds)}`' for duration_seconds, mtype in caps]
        embed = discord.Embed(
            title=f'Active Caps for {channel_obj.mention}',
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
        if interaction.guild is None:
            return await interaction.response.send_message(content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, target)
        found_aliases = False
        lines = []
        highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content='\U0001F6AB Only owners, developers or administrators can list all aliases across the server.')
            aliases = await Alias.fetch_command_aliases_by_guild(interaction.guild)
            if not aliases:
                return await interaction.response.send_message(content=f'\U0001F6AB No aliases found.')
            lines.extend(Alias.format_aliases(aliases))
            found_aliases = len(lines) > 0
            if not found_aliases:
                return await interaction.response.send_message(content='\U0001F6AB No aliases found in this server.')
            pages = []
            chunk_size = 18
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(title='All Aliases in Server', description='\n'.join(chunk), color=discord.Color.blue())
                pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        else:
            if channel_obj is None or channel_obj.type != discord.ChannelType.voice:
                return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
            aliases = await Alias.fetch_command_aliases_by_channel(channel_obj)
            if not aliases:
                return await interaction.response.send_message(content=f'\U0001F6AB No aliases found.')
            lines.extend(Alias.format_aliases(aliases))
            found_aliases = len(lines) > 0
            if not found_aliases:
                return await interaction.response.send_message(content=f'\U0001F6AB No aliases found for the requested target: {target if target else channel_obj.mention}.')
        await interaction.response.send_message(embed=embed)
        
    # DONE
    @commands.command(name='cmds', help='List command aliases routed to a specific channel, temp room, or all channels if "all" is provided.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_commands_text_command(self, ctx: commands.Context, target: Optional[str] = commands.parameter(default=None, description='"all", channel name/ID/mention, or temp room name')) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        found_aliases = False
        lines = []
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content='\U0001F6AB Only owners, developers and administrators can list all aliases across the server.')
            aliases = await Alias.fetch_command_aliases_by_guild(ctx.guild)
            lines.extend(Alias.format_aliases(aliases))
            found_aliases = len(lines) > 0
            if not found_aliases:
                return await self.handler.send_message(ctx, content='\U0001F6AB No aliases found in this server.')
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
            aliases = await Alias.fetch_command_aliases_by_channel(channel_obj)
            if not aliases:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No aliases found.')
            lines.extend(Alias.format_aliases(aliases))
            found_aliases = len(lines) > 0
            if not found_aliases:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No aliases found for the requested target: {target if target else channel_obj.mention}.')
        embed = discord.Embed(title=f'Aliases for {channel_obj.mention}', description='\n'.join(lines), color=discord.Color.blue())
        await self.handler.send_message(ctx, embed=embed)
 
    # DONE
    @app_commands.command(name='del', description='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @app_commands.describe(message_id='Message snowflake ID', channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def delete_message_app_command(
        self,
        interaction: discord.Interaction,
        message_id: Optional[str],
        channel: Optional[str]
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        try:
            msg = await channel_obj.fetch_message(int(message_id))
        except:
            logger.warning('No message with this ID exists.')
        if not msg:
            return await interaction.response.send_message(content=f'\U0001F6AB No message with ID `{message_id}` found in `{channel_obj.name}`.')
        try:
            await msg.delete()
            return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} Message `{message_id}` deleted successfully.')
        except discord.Forbidden:
            logger.warning('Missing permissions to delete the message.')
        return await interaction.response.send_message(content='\U0001F6AB Failed to delete the message.')

    # DONE
    @commands.command(name='del', help='Delete a message by ID (only if you are coordinator/moderator of that temp room).')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def delete_message_text_command(
        self,
        ctx: commands.Context,
        message_id: Optional[int] = commands.parameter(default=None, description='Message snowflake'),
        *,
        channel: Optional[str] = commands.parameter(default=None, description='Channel or snowflake')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        try:
            msg = await channel_obj.fetch_message(message_id)
        except:
            logger.warning('No message with this ID exists.')
        if not msg:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB No message with ID `{message_id}` found in `{channel_obj.name}`.')
        try:
            await msg.delete()
            return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} Message `{message_id}` deleted successfully.')
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
            if not rows: return await interaction.response.send_message(content='\U0001F6AB No users are flagged in any voice channels.')
            channel_map = defaultdict(list)
            for row in rows: channel_map[row['channel_id']].append(row['discord_snowflake'])
            pages = []
            for ch_id, user_ids in sorted(channel_map.items()):
                ch = interaction.guild.get_channel(ch_id)
                ch_name = ch.mention if ch else f'Unknown Channel ({ch_id})'
                embed = discord.Embed(
                    title=f'\U0001F6A9 Flag records for {ch_name}',
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
                return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} is not flagged in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
            lines = []
            for r in rows:
                ch = interaction.guild.get_channel(r['channel_id'])
                ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                reason = r['reason'] or "No reason given"
                lines.append(f'• {ch_name} — {reason}')
            embed=discord.Embed(
                title=f'\U0001F6A9 Flag records for {member_obj.display_name}',
                description='\n'.join(lines),
                color=discord.Color.orange()
            )
            return await interaction.response.send_message(embed=embed,allowed_mentions=discord.AllowedMentions.none())
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                rows = await conn.fetch('SELECT discord_snowflake FROM active_flags WHERE guild_id=$1 AND channel_id=$2', interaction.guild.id, channel_obj.id)
            if not rows:
                return await interaction.response.send_message(content=f'\U0001F6AB No users are flagged for {channel_obj.mention}.')
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
                        title=f'\U0001F6A9 Flag records for {channel_obj.name}',
                        color=discord.Color.red()
                    )
                    embed.add_field(name=f'{interaction.guild.name}', value='\n'.join(formatted_lines), inline=False)
                    pages.append(embed)
            if not pages:
                return await interaction.response.send_message(content=f'\U0001F6AB No flagged users currently in {interaction.guild.name}.')
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
                    title=f'\U0001F6A9 Flag records for {ch_name}',
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
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not flagged in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in rows:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    reason = r['reason'] or "No reason given"
                    lines.append(f'• {ch_name} — {reason}')
                embed = discord.Embed(
                    title=f'\U0001F6A9 Flag records for {member_obj.display_name}',
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
                        title=f'\U0001F6A9 Flag records for {channel_obj.name}',
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
                    return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} is not cowed in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in rows:
                    ch = interaction.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    created_at = discord.utils.format_dt(r['created_at'],style='R') if r['created_at'] else ''
                    lines.append(f'• {ch_name} — {created_at}')
                embed=discord.Embed(
                    title=f'🐮 {member_obj.display_name}',
                    description='\n'.join(lines),
                    color=discord.Color.green()
                )
                return await interaction.response.send_message(embed=embed, allowed_mentions=discord.AllowedMentions.all())
            elif channel_obj:
                if channel_obj.type != discord.ChannelType.voice:
                    return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
                rows = await conn.fetch('''SELECT discord_snowflake, created_at FROM active_cows WHERE guild_id=$1 AND channel_id=$2''', interaction.guild.id, channel_obj.id)
                if not rows:
                    return await interaction.response.send_message(content=f'\U0001F6AB No users are cowed in {channel_obj.mention}.')
                lines = []
                for row in rows:
                    uid = row['discord_snowflake']
                    m = interaction.guild.get_member(uid)
                    if not m:
                        continue
                    created_at = discord.utils.format_dt(row['created_at'],style='R') if row['created_at'] else ''
                    lines.append(f'• {m.display_name} — <@{uid}> — {created_at}')
                if not lines:
                    return await interaction.response.send_message(content=f'\U0001F6AB No new vegans currently in {interaction.guild.name}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk=lines[i:i+chunk_size]
                    embed=discord.Embed(
                        title=f'🐮 Vegan records for {channel_obj.name}',
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
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not cowed in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in rows:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'`{r["channel_id"]}`'
                    created_at = discord.utils.format_dt(r['created_at'], style='R') if r['created_at'] else ''
                    lines.append(f'• {ch_name} — {created_at}')
                embed = discord.Embed(
                    title=f'🐮 {member_obj.display_name}',
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
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are cowed in {channel_obj.mention}.')
                lines = []
                for row in rows:
                    uid = row['discord_snowflake']
                    m = ctx.guild.get_member(uid)
                    if not m:
                        continue
                    created_at = discord.utils.format_dt(row['created_at'], style='R') if row['created_at'] else ''
                    lines.append(f'• {m.display_name} — <@{uid}> — {created_at}')
                if not lines:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No new vegans currently in {ctx.guild.name}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(lines), chunk_size):
                    chunk = lines[i:i + chunk_size]
                    embed = discord.Embed(
                        title=f'🐮 Vegan records for {channel_obj.name}',
                        description='\n'.join(chunk),
                        color=discord.Color.green()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()

    # DONE
    @app_commands.command(name='mutes', description='Lists mute statistics.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def list_mutes_app_command(
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
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await interaction.response.send_message(content=f'\U0001F6AB You do not have permission to use this command (`mutes`) in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''SELECT discord_snowflake, channel_id, room_name, expires_at, COALESCE(reason,'No reason provided') AS reason FROM active_voice_mutes WHERE guild_id=$1 AND target='user' ORDER BY channel_id, room_name, discord_snowflake''', interaction.guild.id)
            if not records:
                return await interaction.response.send_message(content=f'\U0001F6AB No muted users currently in {interaction.guild.name}.')
            grouped = defaultdict(list)
            for r in records:
                grouped[r['channel_id']].append(r)
            pages = []
            for channel_id, user_entries in sorted(grouped.items()):
                channel = interaction.guild.get_channel(channel_id)
                channel_name = channel.mention if channel else f'Unknown Channel ({channel_id})'
                chunk_size = 18
                for i in range(0, len(user_entries), chunk_size):
                    embed = discord.Embed(title=f'\U0001F507 Mutes records for {channel_name}', color=discord.Color.orange())
                    for r in user_entries[i:i + chunk_size]:
                        user_id = r['discord_snowflake']
                        member = interaction.guild.get_member(user_id)
                        name = member.display_name if member else f'User ID {user_id}'
                        mention = member.mention if member else f'`{user_id}`'
                        duration_str = Duration.output_display_from_datetime(r['expires_at'])
                        embed.add_field(name=name, value=f'{mention}\nReason: {r["reason"]}\nDuration: {duration_str}', inline=False)
                    pages.append(embed)
            paginator = AppPaginator(self.bot, interaction, pages)
            return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''SELECT guild_id, channel_id, expires_at, reason FROM active_voice_mutes WHERE discord_snowflake=$1 AND guild_id=$2 AND target='user' AND room_name=$3''', member_obj.id, interaction.guild.id, channel_obj.name)
            records = [r for r in records if interaction.guild.get_channel(r['channel_id'])]
            if not records:
                return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} is not muted in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
            lines = []
            for r in records:
                ch = interaction.guild.get_channel(r['channel_id'])
                ch_mention = ch.mention if ch else f'`{r["channel_id"]}`'
                lines.append(f'• {ch_mention} — {r["reason"]} — {Duration.output_display_from_datetime(r["expires_at"])}')
            embed = discord.Embed(
                title=f'\U0001F507 Mute records for {member_obj.display_name}',
                description='\n'.join(lines),
                color=discord.Color.orange()
            )
            return await interaction.response.send_message(embed=embed, allowed_mentions=discord.AllowedMentions.none())
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''SELECT guild_id, channel_id, expires_at, reason, discord_snowflake FROM active_voice_mutes WHERE channel_id=$1 AND guild_id=$2 AND target='user' AND room_name=$3''', channel_obj.id, interaction.guild.id, channel_obj.name)
            if not records:
                return await interaction.response.send_message(content=f'\U0001F6AB No users are currently muted in {channel_obj.mention}.')
            lines = []
            for r in records:
                uid = r['discord_snowflake']
                member_obj = interaction.guild.get_member(uid)
                if not member_obj:
                    continue
                lines.append(f'• {member_obj.display_name} — <@{uid}> — {Duration.output_display_from_datetime(r["expires_at"])}')
            if not lines:
                return await interaction.response.send_message(content=f'\U0001F6AB No muted users currently in {interaction.name}.')
            chunk_size = 18
            pages = []
            for i in range(0, len(lines), chunk_size):
                chunk = lines[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'\U0001F507 Mute records for {channel_obj.name}',
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
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, target)
        if member_obj:
            target = None
        channel_obj = await self.channel_service.resolve_channel(ctx, target)
        highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
        if target and target.lower() == 'all':
            if highest_role not in ('Owner', 'Developer', 'Administrator'):
                return await self.handler.send_message(ctx, content=f'\U0001F6AB You do not have permission to use this command (`mutes`) in {channel_obj.mention}.')
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                    SELECT discord_snowflake, channel_id, room_name, expires_at, COALESCE(reason, 'No reason provided') AS reason
                    FROM active_voice_mutes
                    WHERE guild_id = $1
                      AND target = 'user'
                    ORDER BY channel_id, room_name, discord_snowflake
                ''', ctx.guild.id)
            if not records:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No muted users currently in {ctx.guild.name}.')
            grouped = defaultdict(list)
            for record in records:
                grouped[record['channel_id']].append(record)
            pages = []
            for channel_id, user_entries in sorted(grouped.items()):
                channel = ctx.guild.get_channel(channel_id)
                channel_name = channel.mention if channel else f'Unknown Channel ({channel_id})'
                chunk_size = 18
                for i in range(0, len(user_entries), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F507 Mutes records for {channel_name}',
                        color=discord.Color.orange()
                    )
                    for record in user_entries[i:i + chunk_size]:
                        user_id = record['discord_snowflake']
                        member = ctx.guild.get_member(user_id)
                        name = member.display_name if member else f'User ID {user_id}'
                        mention = member.mention if member else f'`{user_id}`'
                        reason = record['reason'] or 'No reason provided'
                        duration_str = Duration.output_display_from_datetime(record['expires_at'])
                        embed.add_field(name=name, value=f'{mention}\nReason: {reason}\nDuration: {duration_str}', inline=False)
                    pages.append(embed)
            paginator = Paginator(self.bot, ctx, pages)
            return await paginator.start()
        if member_obj:
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                    SELECT guild_id, channel_id, expires_at, reason
                    FROM active_voice_mutes
                    WHERE discord_snowflake = $1
                      AND guild_id = $2
                      AND target = 'user'
                      AND room_name = $3
                ''', member_obj.id, ctx.guild.id, channel_obj.name)
            records = [r for r in records if ctx.guild.get_channel(r['channel_id'])]
            if not records:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not muted in any voice channels.', allowed_mentions=discord.AllowedMentions.none())
            description_lines = []
            for record in records:
                channel_obj = ctx.guild.get_channel(record['channel_id'])
                channel_mention = channel_obj.mention if channel_obj else f'`{record['channel_id']}`'
                reason = record['reason']
                duration_str = Duration.output_display_from_datetime(record['expires_at'])
                description_lines.append(f'• {channel_mention} — {reason} — {duration_str}')
            embed = discord.Embed(title=f'\U0001F507 Mute records for {member_obj.display_name}', description='\n'.join(description_lines), color=discord.Color.orange())
            return await self.handler.send_message(ctx, embed=embed, allowed_mentions=discord.AllowedMentions.none())
        elif channel_obj:
            if channel_obj.type != discord.ChannelType.voice:
                return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
            async with self.bot.db_pool.acquire() as conn:
                records = await conn.fetch('''
                    SELECT guild_id, channel_id, expires_at, reason, discord_snowflake
                    FROM active_voice_mutes
                    WHERE channel_id = $1
                      AND guild_id = $2
                      AND target = 'user'
                      AND room_name = $3
                ''', channel_obj.id, ctx.guild.id, channel_obj.name)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are currently muted in {channel_obj.mention}.')
                description_lines = []
                for record in records:
                    uid = record['discord_snowflake']
                    member_obj = ctx.guild.get_member(uid)
                    if not member_obj:
                        continue
                    name = member_obj.display_name
                    duration_str = Duration.output_display_from_datetime(record['expires_at'])
                    description_lines.append(f'• {name} — <@{uid}> — {duration_str}')
                if not description_lines:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No muted users currently in {channel_obj.mention}.')
                chunk_size = 18
                pages = []
                for i in range(0, len(description_lines), chunk_size):
                    chunk = description_lines[i:i + chunk_size]
                    embed = discord.Embed(title=f'\U0001F507 Mute records for {channel_obj.name}', color=discord.Color.orange())
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
        member: Optional[str] = None,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        member_obj = await self.member_service.resolve_member(interaction, member)
        if not member_obj:
            return await interaction.response.send_message(content=f'\U0001F6AB Invalid member for {member}.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''
                SELECT channel_id, room_name
                FROM active_stages
                WHERE guild_id=$1
            ''', interaction.guild.id)
            if not records:
                return await interaction.response.send_message(content='\U0001F6AB No active stage found.')
            stage = await conn.fetchrow('''
                SELECT initiator_id
                FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', interaction.guild.id, channel_obj.id, channel_obj.name)
            if not stage:
                return await interaction.response.send_message(content='\U0001F6AB No active stage found.')
            if member_obj.id == stage['initiator_id']:
                return await interaction.response.send_message(content='\U0001F6AB You cannot mute the stage initiator.')
            highest_role = await is_owner_developer_administrator_coordinator_moderator(interaction)
            is_coordinator = await conn.fetchval('''
                SELECT 1
                FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
            ''', interaction.guild.id, channel_obj.id, channel_obj.name, interaction.user.id)
            if not is_coordinator and highest_role != 'Everyone':
                return await interaction.response.send_message(content='\U0001F6AB Only stage coordinators or above can use this command.')
            try:
                await member_obj.edit(mute=not member_obj.voice.mute)
                return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been {"muted" if member_obj.voice.mute else "unmuted"}.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Failed to toggle mute: {e}')
            return await interaction.response.send_message(content=f'\U0001F6AB Failed to toggle mute for {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                
    # DONE
    @commands.command(name='mstage', help='Mute/unmute a member in the active stage.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def stage_mute_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: Optional[str] = commands.parameter(default=None, description="Tag a channel or include it's snowflake ID")
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from target: `{member}`.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        async with self.bot.db_pool.acquire() as conn:
            records = await conn.fetch('''
                SELECT channel_id, room_name
                FROM active_stages
                WHERE guild_id=$1
            ''', ctx.guild.id)
            if not records:
                return await self.handler.send_message(ctx, content='\U0001F6AB No active stage found.')
            stage = await conn.fetchrow('''
                SELECT initiator_id
                FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3
            ''', ctx.guild.id, channel_obj.id, channel_obj.name)
            if not stage:
                return await self.handler.send_message(ctx, content='\U0001F6AB No active stage found.')
            if member_obj.id == stage['initiator_id']:
                return await self.handler.send_message(ctx, content='\U0001F6AB You cannot mute the stage initiator.')
            highest_role = await is_owner_developer_administrator_coordinator_moderator(ctx)
            is_coordinator = await conn.fetchval('''
                SELECT 1
                FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
            ''', ctx.guild.id, channel_obj.id, channel_obj.name, ctx.author.id)
            if not is_coordinator and highest_role != 'Everyone':
                return await self.handler.send_message(ctx, content='\U0001F6AB Only stage coordinators can use this command.')
            try:
                await member_obj.edit(mute=not member_obj.voice.mute)
                return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been {"muted" if member_obj.voice.mute else "unmuted"}.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Failed to toggle mute: {e}')
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Failed to toggle mute for {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
    
    # DONE
    @app_commands.command(name='pstage', description='Promote/demote a member as stage coordinator.')
    @app_commands.describe(member='Tag a member or include their snowflake ID')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def stage_promote_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        member_obj = await self.member_service.resolve_member(interaction, member)
        if not member_obj:
            return await interaction.response.send_message(content=f'\U0001F6AB Could not resolve a valid member from target: `{member}`.')
        channel_obj = await self.channel_service.resolve_channel(interaction, None)
        async with self.bot.db_pool.acquire() as conn:
            stage = await conn.fetchrow('SELECT initiator_id, channel_id, room_name FROM active_stages WHERE guild_id=$1', interaction.guild.id)
            if not stage:
                return await interaction.response.send_message(content='\U000026A0\U0000FE0F No active stage found.')
            if member_obj.id == stage['initiator_id']:
                return await interaction.response.send_message(content='\U0001F6AB Cannot change the initiator role.')
            is_coordinator = await conn.fetchval('''
                SELECT 1 FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
            ''', interaction.guild.id, channel_obj.id, channel_obj.name, member_obj.id)
            try:
                if is_coordinator:
                    await conn.execute('''
                        DELETE FROM stage_coordinators
                        WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
                    ''', interaction.guild.id, channel_obj.id, channel_obj.name, member_obj.id)
                    if member_obj.voice and not member_obj.voice.mute:
                        await member_obj.edit(mute=True)
                    return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been demoted from stage coordinator.', allowed_mentions=discord.AllowedMentions.none())
                else:
                    await conn.execute('''
                        INSERT INTO stage_coordinators (guild_id, channel_id, room_name, discord_snowflake)
                        VALUES ($1,$2,$3,$4) ON CONFLICT DO NOTHING
                    ''', interaction.guild.id, channel_obj.id, channel_obj.name, member_obj.id)
                    await conn.execute('''
                        DELETE FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='room'
                    ''', interaction.guild.id, member_obj.id, channel_obj.id)
                    if member_obj.voice and member_obj.voice.mute:
                        await member_obj.edit(mute=False)
                    return await interaction.response.send_message(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been promoted to stage coordinator.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Failed to toggle promotion: {e}')
            return await interaction.response.send_message(content=f'\U0001F6AB Failed to toggle promotion for {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
    
    # DONE
    @commands.command(name='pstage', help='Promote/demote a member as stage coordinator.')
    @is_owner_developer_administrator_coordinator_moderator_predicator()
    async def stage_promote_text_command(
        self,
        ctx: commands.Context,
        member: Optional[str] = commands.parameter(default=None, description='Tag a member or include their snowflake ID')
    ) -> None:
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(ctx, member)
        if not member_obj:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB Could not resolve a valid member from target: `{member}`.')
        channel_obj = await self.channel_service.resolve_channel(ctx, None)
        async with self.bot.db_pool.acquire() as conn:
            stage = await conn.fetchrow('SELECT initiator_id, channel_id, room_name FROM active_stages WHERE guild_id=$1', ctx.guild.id)
            if not stage:
                return await self.handler.send_message(ctx, content='\U0001F6AB No active stage found.')
            if member_obj.id == stage['initiator_id']:
                return await self.handler.send_message(ctx, content='\U0001F6AB Cannot change the initiator role.')
            is_coordinator = await conn.fetchval('''
                SELECT 1 FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
            ''', ctx.guild.id, channel_obj.id, channel_obj.name, member_obj.id)
            try:
                if is_coordinator:
                    await conn.execute('''
                        DELETE FROM stage_coordinators
                        WHERE guild_id=$1 AND channel_id=$2 AND room_name=$3 AND discord_snowflake=$4
                    ''', ctx.guild.id, channel_obj.id, channel_obj.name, member_obj.id)
                    if member_obj.voice and not member_obj.voice.mute:
                        await member_obj.edit(mute=True)
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been demoted from stage coordinator.', allowed_mentions=discord.AllowedMentions.none())
                    
                else:
                    await conn.execute('''
                        INSERT INTO stage_coordinators (guild_id, channel_id, room_name, discord_snowflake)
                        VALUES ($1,$2,$3,$4) ON CONFLICT DO NOTHING
                    ''', ctx.guild.id, channel_obj.id, channel_obj.name, member_obj.id)
                    await conn.execute('''
                        DELETE FROM active_voice_mutes
                        WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND target='room'
                    ''', ctx.guild.id, member_obj.id, channel_obj.id)
                    if member_obj.voice and member_obj.voice.mute:
                        await member_obj.edit(mute=False)
                    return await self.handler.send_message(ctx, content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been promoted to stage coordinator.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.warning(f'Failed to toggle promotion: {e}')
                return await self.handler.send_message(ctx, content=f'\U0001F6AB Failed to toggle promotion for {member_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

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
                stages = await conn.fetch('''
                    SELECT s.channel_id, s.room_name, s.initiator_id, s.expires_at, COUNT(v.discord_snowflake) AS active_mutes
                    FROM active_stages s LEFT JOIN active_voice_mutes v ON s.guild_id=v.guild_id AND s.channel_id=v.channel_id AND v.target='room'
                    WHERE s.guild_id=$1 GROUP BY s.channel_id, s.room_name, s.initiator_id, s.expires_at ORDER BY s.channel_id
                ''', interaction.guild.id)
                if not stages:
                    return await interaction.response.send_message(content=f'\U0001F6AB No active stages in {interaction.guild.name}.')
                pages, chunk_size = [], 8
                for i in range(0, len(stages), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F399 Active Stages in {interaction.guild.name}',
                        color=discord.Color.purple()
                    )
                    for s in stages[i:i+chunk_size]:
                        ch = interaction.guild.get_channel(s['channel_id'])
                        ch_name = ch.mention if ch else f'Unknown Channel ({s["channel_id"]})'
                        embed.add_field(name=ch_name, value=f'Active stage mutes: {s["active_mutes"]}', inline=False)
                    pages.append(embed)
                paginator = AppPaginator(self.bot, interaction, pages)
                return await paginator.start()
            stage = await conn.fetchrow('''
                SELECT initiator_id, expires_at FROM active_stages WHERE guild_id=$1 AND channel_id=$2 AND room_name = $3
            ''', interaction.guild.id, channel_obj.id, channel_obj.name)
            if not stage:
                return await interaction.response.send_message(content=f'\U0001F6AB No active stage in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            mutes = await conn.fetch('''
                SELECT discord_snowflake, expires_at, reason FROM active_voice_mutes WHERE guild_id=$1 AND channel_id=$2 AND target='room' AND room_name = $3
            ''', interaction.guild.id, channel_obj.id, channel_obj.name)
            coordinators = await conn.fetch('''
                SELECT discord_snowflake FROM stage_coordinators WHERE guild_id=$1 AND channel_id=$2
            ''', interaction.guild.id, channel_obj.id)
            coordinator_mentions = [member.mention for c in coordinators if (member:=interaction.guild.get_member(c['discord_snowflake']))]
            coordinator_str = ', '.join(coordinator_mentions) if coordinator_mentions else 'No coordinators'
            initiator = interaction.guild.get_member(stage['initiator_id'])
            initiator_name = initiator.mention if initiator else f'`{stage["initiator_id"]}`'
            expires = Duration.output_display_from_datetime(stage['expires_at']) if stage['expires_at'] else 'No expiration'
            lines = []
            for m in mutes:
                user = interaction.guild.get_member(m['discord_snowflake'])
                duration_str = Duration.output_display_from_datetime(m['expires_at']) if m['expires_at'] else 'No expiration'
                reason = m['reason'] or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = f'Initiator: {initiator_name}\nStage Expires: {expires}\nCoordinators: {coordinator_str}\nActive stage mutes: {len(lines)}\n\n'+'\n'.join(lines)
            pages, chunk_size = [], 18
            for i in range(0, len(description.splitlines()), chunk_size):
                embed=discord.Embed(
                    title=f'\U0001F399 Stage info for {channel_obj.mention}',
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
                stages = await conn.fetch('''
                    SELECT s.channel_id, s.initiator_id, s.expires_at, COUNT(v.discord_snowflake) AS active_mutes
                    FROM active_stages s
                    LEFT JOIN active_voice_mutes v
                        ON s.guild_id=v.guild_id AND s.channel_id=v.channel_id AND v.target='room'
                    WHERE s.guild_id=$1 AND s.room_name = $2
                    GROUP BY s.channel_id, s.initiator_id, s.expires_at
                    ORDER BY s.channel_id
                ''', ctx.guild.id, channel_obj.name)
                if not stages:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No active stages in {ctx.guild.name}.')
                pages, chunk_size = [], 8
                for i in range(0, len(stages), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F399 Active Stages in {ctx.guild.name}',
                        color=discord.Color.purple()
                    )
                    for s in stages[i:i+chunk_size]:
                        ch = ctx.guild.get_channel(s['channel_id'])
                        ch_name = ch.mention if ch else f'Unknown Channel ({s["channel_id"]})'
                        embed.add_field(name=ch_name, value=f'Active stage mutes: {s["active_mutes"]}', inline=False)
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            stage = await conn.fetchrow('''
                SELECT initiator_id, expires_at
                FROM active_stages
                WHERE guild_id=$1 AND channel_id=$2 AND room_name = $3
            ''', ctx.guild.id, channel_obj.id, channel_obj.name)
            if not stage:
                return await self.handler.send_message(ctx, content=f'\U0001F6AB No active stage in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            mutes = await conn.fetch('''
                SELECT discord_snowflake, expires_at, reason
                FROM active_voice_mutes
                WHERE guild_id=$1 AND channel_id=$2 AND target='room' and room_name = $3
            ''', ctx.guild.id, channel_obj.id, channel_obj.name)
            coordinators = await conn.fetch('''
                SELECT discord_snowflake
                FROM stage_coordinators
                WHERE guild_id=$1 AND channel_id=$2
            ''', ctx.guild.id, channel_obj.id)
            coordinator_mentions = []
            for c in coordinators:
                member = ctx.guild.get_member(c['discord_snowflake'])
                if member:
                    coordinator_mentions.append(member.mention)
            coordinator_str = ', '.join(coordinator_mentions) if coordinator_mentions else 'No coordinators'
            initiator = ctx.guild.get_member(stage['initiator_id'])
            initiator_name = initiator.mention if initiator else f'`{stage["initiator_id"]}`'
            expires = Duration.output_display_from_datetime(stage['expires_at']) if stage['expires_at'] else 'No expiration'
            lines = []
            for m in mutes:
                user = ctx.guild.get_member(m['discord_snowflake'])
                duration_str = Duration.output_display_from_datetime(m['expires_at']) if m['expires_at'] else 'No expiration'
                reason = m['reason'] or 'No reason provided'
                lines.append(f'• {user.mention} — {reason} — {duration_str}')
            description = (
                f'Initiator: {initiator_name}\n'
                f'Stage Expires: {expires}\n'
                f'Coordinators: {coordinator_str}\n'
                f'Active stage mutes: {len(lines)}\n\n' + '\n'.join(lines)
            )
            pages, chunk_size = [], 18
            for i in range(0, len(description.splitlines()), chunk_size):
                embed = discord.Embed(
                    title=f'\U0001F399 Stage info for {channel_obj.mention}',
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
                records = await conn.fetch('''SELECT discord_snowflake, channel_id, room_name, reason, expires_at FROM active_text_mutes WHERE guild_id = $1 ORDER BY room_name, channel_id, discord_snowflake''', interaction.guild.id)
                if not records:
                    return await interaction.response.send_message(content=f'\U0001F6AB No users are currently text-muted in {interaction.guild.name}.')
                grouped = defaultdict(list)
                for r in records:
                    grouped[r['channel_id']].append(r)
                pages, chunk_size = [], 18
                for ch_id, entries in sorted(grouped.items()):
                    ch = interaction.guild.get_channel(ch_id)
                    ch_name = ch.mention
                    for i in range(0, len(entries), chunk_size):
                        embed= discord.Embed(
                            title=f'\U0001F4DA Text-mute records for {ch_name}',
                            color=discord.Color.orange()
                        )
                        for e in entries[i:i+chunk_size]:
                            user = interaction.guild.get_member(e['discord_snowflake'])
                            mention = user.name if user else f'`{e["discord_snowflake"]}`'
                            reason = e['reason'] or 'No reason provided'
                            duration_str = Duration.output_display_from_datetime(e['expires_at'])
                            embed.add_field(name=mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
                        pages.append(embed)
                paginator = AppPaginator(self.bot, interaction, pages)
                return await paginator.start()
            elif member_obj:
                records = await conn.fetch('''SELECT channel_id, room_name, reason, expires_at FROM active_text_mutes WHERE discord_snowflake = $1 AND guild_id = $2 AND room_name = $3''', member_obj.id, interaction.guild.id, channel_obj.name)
                if not records: return await interaction.response.send_message(content=f'\U0001F6AB {member_obj.mention} is not text-muted in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    ch = interaction.guild.get_channel(r['channel_id'])
                    duration_str = Duration.output_display_from_datetime(r['expires_at'])
                    lines.append(f'• {ch.mention} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F4DA Text-mute records for {member_obj.display_name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = AppPaginator(self.bot, interaction, pages)
                return await paginator.start()
            elif channel_obj:
                records = await conn.fetch('''SELECT discord_snowflake, room_name, reason, expires_at FROM active_text_mutes WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3''', channel_obj.id, interaction.guild.id, channel_obj.name)
                if not records:
                    return await interaction.response.send_message(content=f'\U0001F6AB No users are currently text-muted in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    user = interaction.guild.get_member(r['discord_snowflake'])
                    if not user:
                        continue
                    duration_str = Duration.output_display_from_datetime(r['expires_at'])
                    lines.append(f'• {user.mention} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F4DA Text-mute records for {channel_obj.name}',
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
                records = await conn.fetch('''
                    SELECT discord_snowflake, channel_id, room_name, reason, expires_at
                    FROM active_text_mutes
                    WHERE guild_id = $1
                    ORDER BY room_name, channel_id, discord_snowflake
                ''', ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are currently text-muted in {ctx.guild.name}.')
                grouped = defaultdict(list)
                for r in records:
                    key = (r['channel_id'], r['room_name'])
                    grouped[key].append(r)
                pages, chunk_size = [], 18
                for (ch_id, rname), entries in sorted(grouped.items()):
                    ch = ctx.guild.get_channel(ch_id)
                    ch_name = ch.mention if ch else f'Channel ID `{ch_id}`'
                    if rname:
                        ch_name = f"{ch_name}"
                    for i in range(0, len(entries), chunk_size):
                        embed = discord.Embed(
                            title=f'\U0001F4DA Text-mute records for {ch_name}',
                            color=discord.Color.orange()
                        )
                        for e in entries[i:i + chunk_size]:
                            user = ctx.guild.get_member(e['discord_snowflake'])
                            mention = user.name if user else f'<@{e["discord_snowflake"]}>'
                            reason = e['reason'] or 'No reason provided'
                            duration_str = Duration.output_display_from_datetime(e['expires_at'])
                            embed.add_field(name=mention, value=f'Reason: {reason}\nDuration: {duration_str}', inline=False)
                        pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif member_obj:
                records = await conn.fetch('''
                    SELECT channel_id, room_name, reason, expires_at
                    FROM active_text_mutes
                    WHERE discord_snowflake = $1 AND guild_id = $2
                ''', member_obj.id, ctx.guild.id)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB {member_obj.mention} is not text-muted in any channels.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    ch = ctx.guild.get_channel(r['channel_id'])
                    ch_name = ch.mention if ch else f'Channel ID `{r["channel_id"]}`'
                    if r['room_name']:
                        ch_name = f"{ch_name}"
                    duration_str = Duration.output_display_from_datetime(r['expires_at'])
                    lines.append(f'• {ch_name} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F4DA Text-mute records for {member_obj.display_name}',
                        description='\n'.join(lines[i:i+chunk_size]),
                        color=discord.Color.orange()
                    )
                    pages.append(embed)
                paginator = Paginator(self.bot, ctx, pages)
                return await paginator.start()
            elif channel_obj:
                records = await conn.fetch('''
                    SELECT discord_snowflake, reason, expires_at
                    FROM active_text_mutes
                    WHERE channel_id = $1 AND guild_id = $2 AND room_name = $3
                ''', channel_obj.id, ctx.guild.id, channel_obj.name)
                if not records:
                    return await self.handler.send_message(ctx, content=f'\U0001F6AB No users are currently text-muted in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
                lines = []
                for r in records:
                    user = ctx.guild.get_member(r['discord_snowflake'])
                    mention = user.mention if user else f'<@{r["discord_snowflake"]}>'
                    duration_str = Duration.output_display_from_datetime(r['expires_at'])
                    lines.append(f'• {mention} — {r["reason"]} — {duration_str}')
                pages, chunk_size = [], 18
                for i in range(0, len(lines), chunk_size):
                    embed = discord.Embed(
                        title=f'\U0001F4DA Text-mute records for {channel_obj.name}',
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

