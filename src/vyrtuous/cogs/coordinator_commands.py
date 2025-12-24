''' coordinator_commands.py A discord.py cog containing coordinator commands for the Vyrtuous bot.

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
from vyrtuous.service.discord_message_service import DiscordMessageService
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.check_service import *
from vyrtuous.utils.administrator import Administrator
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.moderator import Moderator

class CoordinatorCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.handler = DiscordMessageService(self.bot, self.bot.db_pool)
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.member_service = MemberService()
        
    # DONE
    @app_commands.command(name='mod', description="Grants/revokes a user's permission to `Moderator` for a specific channel.")
    @app_commands.describe(member='Tag a member or include their snowflake ID', channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_coordinator_predicator()
    async def create_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: Optional[str] = None,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(interaction, member)
        if not member_obj:
            return await interaction.response.send_message(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.id == interaction.guild.me.id:
            return await interaction.response.send_message(content='\U0001F6AB You cannot make the bot a moderator.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        channel_related_role = await is_owner_developer_administrator_coordinator_via_channel_member(channel_obj, interaction.user)
        if channel_related_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
            return await interaction.response.send_message(content=f'\U0001F6AB You are not permitted to grant/revoke moderator status for {channel_obj.mention}..')
        success = await has_equal_or_higher_role(interaction, member_obj, channel_obj)
        if not success:
            return await interaction.response.send_message(content=f'\U0001F6AB You are not allowed to add/remove {member_obj.mention} as a moderator because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            action = None
            moderator_channel_ids = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            if moderator_channel_ids and channel_obj.id in moderator_channel_ids:
                await Moderator.delete_by_channel_and_member(channel_snowflake=channel_obj.id, member_snowflake=member_obj.id)
                action = 'revoked'
            else:
                moderator = Moderator(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
                await moderator.create()
                action = 'granted'
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_moderator', member_obj.id, interaction.user.id, interaction.guild.id, channel_obj.id, f'Moderator access {action}')
        return await interaction.response.send_message(content=f"{self.emoji.get_random_emoji()} {member_obj.mention}'s moderator access has been {action} in {channel_obj.mention}")
              
    # DONE
    @commands.command(name='mod', help='Grants/revokes a user to `Moderator` for a specific channel.')
    @is_owner_developer_administrator_coordinator_predicator()
    async def create_moderator_text_command(
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
        if member_obj.id == ctx.guild.me.id:
            return await self.handler.send_message(ctx, content='\U0001F6AB You cannot make the bot a moderator.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        channel_related_role = await is_owner_developer_administrator_coordinator_via_channel_member(channel_obj, ctx.author)
        if channel_related_role not in ('Owner', 'Developer', 'Administrator', 'Coordinator'):
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not permitted to grant/revoke moderator status for {channel_obj.mention}.')
        success = await has_equal_or_higher_role(ctx.message, member_obj, channel_obj)
        if not success:
            return await self.handler.send_message(ctx, content=f'\U0001F6AB You are not allowed to add/remove {member_obj.mention} as a moderator because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            action = None
            moderator_channel_ids = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            if moderator_channel_ids and channel_obj.id in moderator_channel_ids:
                await Moderator.delete_by_channel_and_member(channel_snowflake=channel_obj.id, member_snowflake=member_obj.id)
                action = 'revoked'
            else:
                moderator = Moderator(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
                await moderator.create()
                action = 'granted'
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'toggled_moderator', member_obj.id, ctx.author.id, ctx.guild.id, channel_obj.id, f'Moderator access {action}')
        return await self.handler.send_message(ctx, content=f"{self.emoji.get_random_emoji()} {member_obj.mention}'s moderator access has been {action} in {channel_obj.mention}.")
 
    # DONE
    @app_commands.command(name='rmute', description='Mutes all members in a VC (except yourself).')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_coordinator_predicator()
    async def room_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='\U0001F6AB This command can only be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        muted_members, skipped_members, failed_members = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel_obj.members:
                if member.id == interaction.user.id:
                    continue
                try:
                    voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target="user")
                    if voice_mute:
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=True)
                    voice_mute = VoiceMute(channel_snowflake=channel_obj.id, expires_at=None, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target="user")
                    voice_mute.create()
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'mute', member.id, interaction.user.id, interaction.guild.id, channel_obj.id, 'Muted via room_mute')
                    muted_members.append(member)
                except Exception as e:
                    logger.warning(f'Mute failed for {member.mention}: {e}')
                    failed_members.append(member)
        summary = f'{self.emoji.get_random_emoji()} Muted {len(muted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to mute {len(failed_members)}.'
        return await interaction.response.send_message(content=summary)

    # DONE
    @commands.command(name='rmute', help='Mutes all members in a VC (except yourself).')
    @is_owner_developer_administrator_coordinator_predicator()
    async def room_mute_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        muted_members, skipped_members, failed_members = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel_obj.members:
                if member.id == ctx.author.id:
                    continue
                try:
                    voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target="user")
                    if voice_mute:
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=True)
                    voice_mute = VoiceMute(channel_snowflake=channel_obj.id, expires_at=None, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target="user")
                    voice_mute.create()
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'mute', member.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Muted via room_mute')
                    muted_members.append(member)
                except Exception as e:
                    logger.warning(f'Mute failed for {member.mention}: {e}')
                    failed_members.append(member)
        summary = f'{self.emoji.get_random_emoji()} Muted {len(muted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to mute {len(failed_members)}.'
        return await self.handler.send_message(ctx, content=summary)

    # DONE
    @app_commands.command(name='xrmute', description='Unmutes all members in a VC (except yourself).')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_coordinator_predicator()
    async def room_unmute_app_command(
        self,
        interaction: discord.Interaction,
        channel: Optional[str] = None
    ):
        if not interaction.guild:
            return await interaction.response.send_message(content='This command must be used in a server.')
        channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await interaction.response.send_message(content='\U0001F6AB Please specify a valid target.')
        unmuted_members, skipped_members, failed_members = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel_obj.members:
                try:
                    voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target="user")
                    if not voice_mute:
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=False)
                    await VoiceMute.delete_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target="user")
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'unmute', member.id, interaction.user.id, interaction.guild.id, channel_obj.id, 'Unmuted via room_unmute')
                    unmuted_members.append(member)
                except Exception as e:
                    logger.warning(f'Unmute failed for {member.mention}: {e}')
                    failed_members.append(member)
        summary = f'{self.emoji.get_random_emoji()} Unmuted {len(unmuted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to unmute {len(failed_members)}.'
        return await interaction.response.send_message(content=summary)
    
    # DONE
    @commands.command(name='xrmute', help='Unmutes all members in a VC (except yourself).')
    @is_owner_developer_administrator_coordinator_predicator()
    async def room_unmute_text_command(
        self,
        ctx: commands.Context,
        channel: Optional[str] = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        if not ctx.guild:
            return await self.handler.send_message(ctx, content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        if channel_obj.type != discord.ChannelType.voice:
            return await self.handler.send_message(ctx, content='\U0001F6AB Please specify a valid target.')
        unmuted_members, skipped_members, failed_members = [], [], []
        async with self.bot.db_pool.acquire() as conn:
            for member in channel_obj.members:
                if member.id == ctx.author.id:
                    continue
                try:
                    voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target="user")
                    if not voice_mute:
                        skipped_members.append(member)
                        continue
                    if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                        await member.edit(mute=False)
                    await VoiceMute.delete_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target="user")
                    await conn.execute('''
                        INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                        VALUES ($1, $2, $3, $4, $5, $6)
                    ''', 'unmute', member.id, ctx.author.id, ctx.guild.id, channel_obj.id, 'Unmuted via room_unmute')
                    unmuted_members.append(member)
                except Exception as e:
                    logger.warning(f'Unmute failed for {member.mention}: {e}')
                    failed_members.append(member)
        summary = f'{self.emoji.get_random_emoji()} Unmuted {len(unmuted_members)} member(s) in {channel_obj.mention}.'
        if skipped_members:
            summary += f'\n\U000026A0\U0000FE0F Skipped {len(skipped_members)}.'
        if failed_members:
            summary += f'\n\U0001F6AB Failed to unmute {len(failed_members)}.'
        return await self.handler.send_message(ctx, content=summary)
    
async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorCommands(bot))
