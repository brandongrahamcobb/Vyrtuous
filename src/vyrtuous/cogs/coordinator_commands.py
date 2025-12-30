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
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.check_service import *
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.moderator import Moderator
from vyrtuous.utils.state import State
from vyrtuous.utils.voice_mute import VoiceMute
from vyrtuous.utils.snowflake import *

class CoordinatorCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.channel_service = ChannelService()
        self.emoji = Emojis()
        self.member_service = MemberService()
        self.message_service = MessageService(self.bot, self.bot.db_pool)
        
    # DONE
    @app_commands.command(name='mod', description="Grants/revokes a user's permission to `Moderator` for a specific channel.")
    @app_commands.describe(member='Tag a member or include their snowflake ID', channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_coordinator_predicator()
    async def create_moderator_app_command(
        self,
        interaction: discord.Interaction,
        member: AppMemberSnowflake,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        action = None
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
            member_obj = await self.member_service.resolve_member(interaction, member)
            check_not_self(interaction, member_snowflake=member_obj.id)
            highest_role = await has_equal_or_higher_role(interaction, channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id, sender_snowflake=interaction.user.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        moderator_channel_ids = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
        if moderator_channel_ids and channel_obj.id in moderator_channel_ids:
            await Moderator.delete_by_channel_and_member(channel_snowflake=channel_obj.id, member_snowflake=member_obj.id)
            action = 'revoked'
        else:
            moderator = Moderator(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member_obj.id)
            await moderator.grant()
            action = 'granted'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Moderator access for {member_obj.mention} has been {action} in {channel_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
              
    # DONE
    @commands.command(name='mod', help='Grants/revokes a user to `Moderator` for a specific channel.')
    @is_owner_developer_administrator_coordinator_predicator()
    async def create_moderator_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(default=None, description='Tag a member or include their snowflake ID'),
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        action = None
        channel_obj = None
        member_obj = None
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
            member_obj = await self.member_service.resolve_member(ctx, member)
            check_not_self(ctx, member_snowflake=member_obj.id)
            highest_role = await has_equal_or_higher_role(ctx, channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id, sender_snowflake=ctx.author.id)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        moderator_channel_ids = await Moderator.fetch_channels_by_guild_and_member(guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
        if moderator_channel_ids and channel_obj.id in moderator_channel_ids:
            await Moderator.delete_by_channel_and_member(channel_snowflake=channel_obj.id, member_snowflake=member_obj.id)
            action = 'revoked'
        else:
            moderator = Moderator(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member_obj.id)
            await moderator.grant()
            action = 'granted'
        try:
            return await state.end(success=f'{self.emoji.get_random_emoji()} Moderator access for {member_obj.mention} has been {action} in {channel_obj.mention}.')
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
        
    # DONE
    @app_commands.command(name='rmute', description='Mutes all members in a VC (except yourself).')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_coordinator_predicator()
    async def room_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        chunk_size = 18
        muted_members, pages, skipped_members, failed_members = [], [], [], []
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        for member in channel_obj.members:
            if member.id == interaction.user.id:
                continue
            voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target='user')
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                try:
                    await member.edit(mute=True)
                except Exception as e:
                    logger.warning(f'Failed to mute {member.id}: {e}.')
                    failed_members.append(member)
            voice_mute = VoiceMute(channel_snowflake=channel_obj.id, expires_at=None, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target='user')
            await voice_mute.create()
            muted_members.append(member)
        muted_chunks = [muted_members[i:i + chunk_size] for i in range(0, len(muted_members), chunk_size)]
        skipped_chunks = [skipped_members[i:i + chunk_size] for i in range(0, len(skipped_members), chunk_size)]
        failed_chunks = [failed_members[i:i + chunk_size] for i in range(0, len(failed_members), chunk_size)]
        member_sets = max(len(muted_chunks), len(skipped_chunks), len(failed_chunks), 1)
        for page in range(member_sets):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Mute Summary',
                color=0xED4245
            )
            if page < len(muted_chunks):
                chunk = muted_chunks[page]
                embed.add_field(
                    name=f'Muted ({len(chunk)}/{len(muted_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            if page < len(skipped_chunks):
                chunk = skipped_chunks[page]
                embed.add_field(
                    name=f'Skipped ({len(chunk)}/{len(skipped_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            if page < len(failed_chunks):
                chunk = failed_chunks[page]
                embed.add_field(
                    name=f'Failed ({len(chunk)}/{len(failed_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            embed.set_footer(text=f'Page {page + 1}/{pages}')
            pages.append(embed)
        try:
             return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @commands.command(name='rmute', help='Mutes all members in a VC (except yourself).')
    @is_owner_developer_administrator_coordinator_predicator()
    async def room_mute_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        chunk_size = 18
        muted_members, pages, skipped_members, failed_members = [], [], [], []
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        for member in channel_obj.members:
            if member.id == ctx.author.id:
                continue
            voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target='user')
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                try:
                    await member.edit(mute=True)
                except Exception as e:
                    logger.warning(f'Failed to mute {member.id}: {e}.')
                    failed_members.append(member)
            voice_mute = VoiceMute(channel_snowflake=channel_obj.id, expires_at=None, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target='user')
            await voice_mute.create()
            muted_members.append(member)
        muted_chunks = [muted_members[i:i + chunk_size] for i in range(0, len(muted_members), chunk_size)]
        skipped_chunks = [skipped_members[i:i + chunk_size] for i in range(0, len(skipped_members), chunk_size)]
        failed_chunks = [failed_members[i:i + chunk_size] for i in range(0, len(failed_members), chunk_size)]
        member_sets = max(len(muted_chunks), len(skipped_chunks), len(failed_chunks), 1)
        for page in range(member_sets):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Mute Summary',
                color=0xED4245
            )
            if page < len(muted_chunks):
                chunk = muted_chunks[page]
                embed.add_field(
                    name=f'Muted ({len(chunk)}/{len(muted_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            if page < len(skipped_chunks):
                chunk = skipped_chunks[page]
                embed.add_field(
                    name=f'Skipped ({len(chunk)}/{len(skipped_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            if page < len(failed_chunks):
                chunk = failed_chunks[page]
                embed.add_field(
                    name=f'Failed ({len(chunk)}/{len(failed_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            embed.set_footer(text=f'Page {page + 1}/{pages}')
            pages.append(embed)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')

    # DONE
    @app_commands.command(name='xrmute', description='Unmutes all members in a VC (except yourself).')
    @app_commands.describe(channel='Tag a channel or include its snowflake ID')
    @is_owner_developer_administrator_coordinator_predicator()
    async def room_unmute_app_command(
        self,
        interaction: discord.Interaction,
        channel: AppChannelSnowflake
    ):
        state = State(interaction)
        channel_obj = None
        chunk_size = 18
        failed_members, pages, skipped_members, unmuted_members = [], [], [], []
        try:
            channel_obj = await self.channel_service.resolve_channel(interaction, channel)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
            except Exception as e:
                return await state.end(error=f'\U0001F3C6 {e}.')
        for member in channel_obj.members:
            voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, member_snowflake=member.id, target='user')
            if not voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                try:
                    await member.edit(mute=False)
                except Exception as e:
                    logger.warning(f'Unmute failed for {member.mention}: {e}.')
                    failed_members.append(member)
            await VoiceMute.delete_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=interaction.guild.id, target='user')         
            unmuted_members.append(member)
        unmuted_chunks = [unmuted_members[i:i + chunk_size] for i in range(0, len(unmuted_members), chunk_size)]
        skipped_chunks = [skipped_members[i:i + chunk_size] for i in range(0, len(skipped_members), chunk_size)]
        failed_chunks = [failed_members[i:i + chunk_size] for i in range(0, len(failed_members), chunk_size)]
        member_sets = max(len(unmuted_chunks), len(skipped_chunks), len(failed_chunks), 1)
        for page in range(member_sets):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Unmute Summary',
                color=0xED4245
            )
            if page < len(unmuted_chunks):
                chunk = unmuted_chunks[page]
                embed.add_field(
                    name=f'Unmuted ({len(chunk)}/{len(unmuted_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            if page < len(skipped_chunks):
                chunk = skipped_chunks[page]
                embed.add_field(
                    name=f'Skipped ({len(chunk)}/{len(skipped_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            if page < len(failed_chunks):
                chunk = failed_chunks[page]
                embed.add_field(
                    name=f'Failed ({len(chunk)}/{len(failed_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            embed.set_footer(text=f'Page {page + 1}/{member_sets}')
            pages.append(embed)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
    # DONE
    @commands.command(name='xrmute', help='Unmutes all members in a VC (except yourself).')
    @is_owner_developer_administrator_coordinator_predicator()
    async def room_unmute_text_command(
        self,
        ctx: commands.Context,
        channel: ChannelSnowflake = commands.parameter(default=None, description='Tag a channel or include its snowflake ID')
    ):
        state = State(ctx)
        channel_obj = None
        chunk_size = 18
        failed_members, pages, skipped_members, unmuted_members = [], [], [], []
        try:
            channel_obj = await self.channel_service.resolve_channel(ctx, channel)
        except Exception as e:
            return await state.end(warning=f'\U000026A0\U0000FE0F {e}.')
        for member in channel_obj.members:
            voice_mute = await VoiceMute.fetch_by_channel_guild_member_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, member_snowflake=member.id, target='user')
            if not voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel and member.voice.channel.id == channel_obj.id:
                try:
                    await member.edit(mute=False)
                except Exception as e:
                    logger.warning(f'Unmute failed for {member.id}: {e}.')
                    failed_members.append(member)
            await VoiceMute.delete_by_channel_guild_and_target(channel_snowflake=channel_obj.id, guild_snowflake=ctx.guild.id, target='user')         
            unmuted_members.append(member)
        unmuted_chunks = [unmuted_members[i:i + chunk_size] for i in range(0, len(unmuted_members), chunk_size)]
        skipped_chunks = [skipped_members[i:i + chunk_size] for i in range(0, len(skipped_members), chunk_size)]
        failed_chunks = [failed_members[i:i + chunk_size] for i in range(0, len(failed_members), chunk_size)]
        member_sets = max(len(unmuted_chunks), len(skipped_chunks), len(failed_chunks), 1)
        for page in range(member_sets):
            embed = discord.Embed(
                title=f'{self.emoji.get_random_emoji()} Unmute Summary',
                color=0xED4245
            )
            if page < len(unmuted_chunks):
                chunk = unmuted_chunks[page]
                embed.add_field(
                    name=f'Unmuted ({len(chunk)}/{len(unmuted_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            if page < len(skipped_chunks):
                chunk = skipped_chunks[page]
                embed.add_field(
                    name=f'Skipped ({len(chunk)}/{len(skipped_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            if page < len(failed_chunks):
                chunk = failed_chunks[page]
                embed.add_field(
                    name=f'Failed ({len(chunk)}/{len(failed_members)})',
                    value='\n'.join(member.mention for member in chunk),
                    inline=False
                )
            embed.set_footer(text=f'Page {page + 1}/{member_sets}')
            pages.append(embed)
        try:
            return await state.end(success=pages)
        except Exception as e:
            return await state.end(error=f'\U0001F3C6 {e}.')
    
async def setup(bot: DiscordBot):
    await bot.add_cog(CoordinatorCommands(bot))
