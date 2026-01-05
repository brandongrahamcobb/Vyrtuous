
''' history.py A utility module for managing and history of messages in the Vyrtuous Discord bot.

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
from datetime import datetime, timezone
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.duration import DurationObject
from vyrtuous.utils.paginator import Paginator 
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.snowflake import *
import discord

class History:
        
    ACTION_TYPES = ['create', 'delete', 'modify']
    ENTRY_TYPES = ['all', 'channel', 'member']

    def __init__(self, channel_snowflake: Optional[int], enabled: Optional[bool], guild_snowflake: Optional[int], snowflakes: list[int|None], entry_type: Optional[str]):
        self.action: Optional[str]
        self.channel_snowflake = channel_snowflake
        self.enabled = enabled
        self.guild_snowflake = guild_snowflake
        self.snowflakes = snowflakes
        self.entry_type = entry_type

    @property
    def action(self):
        self._action

    @action.setter
    def action(self, action: Optional[str]):
        if action not in self.ACTION_TYPES:
            raise ValueError('Invalid action.')
        self._action = action

    @property
    def entry_type(self):
        return self._entry_type

    @entry_type.setter
    def entry_type(self, entry_type: Optional[str]):
        if entry_type not in self.ENTRY_TYPES:
            raise ValueError('Invalid entry type.')
        self._entry_type = entry_type

    @classmethod
    async def send_entry(
        cls,
        alias,
        channel: Optional[discord.VoiceChannel],
        duration: Optional[str],
        highest_role: Optional[str],
        is_channel_scope: bool,
        is_modification: bool,
        member: discord.Member,
        message,
        reason: Optional[str]
    ):
        bot = DiscordBot.get_instance()
        author_snowflake = None
        expires_at = None
        history = await History.fetch_all()
        if message:
            for entry in history:
                channel = bot.get_channel(entry.channel_snowflake)
                pages = cls.build_history_embeds(
                        alias=alias, channel=channel, duration=duration, highest_role=highest_role,
                        is_channel_scope=is_channel_scope, is_modification=is_modification,
                        member=member, message=message,
                        reason=reason
                    )
                paginator = Paginator(bot, channel, pages)
                await paginator.start()
            author_snowflake = message.author.id
        if isinstance(duration, DurationObject):
            expires_at = datetime.now(timezone.utc) + duration.to_timedelta()
        elif duration is not None:
            expires_at = datetime.now(timezone.utc) + DurationObject(duration).to_timedelta()
        else:
            expires_at = datetime.now(timezone.utc) + DurationObject(0).to_timedelta()
        channel_members_voice_count = sum(len(channel.members) for channel in channel.guild.voice_channels)
        guild_members_offline_and_online_member_count = sum(1 for member in channel.guild.members if not member.bot)
        guild_members_online_count = sum(1 for member in channel.guild.members if not member.bot and member.status != discord.Status.offline)
        guild_members_voice_count = sum(len([member for member in channel.members if not member.bot]) for channel in channel.guild.voice_channels)
        await cls.save(action_type=alias.alias_type, channel_members_voice_count=channel_members_voice_count, channel_snowflake=channel.id, executor_member_snowflake=author_snowflake, expires_at=expires_at, guild_members_offline_and_online_member_count=guild_members_offline_and_online_member_count, guild_members_online_count=guild_members_online_count, guild_members_voice_count=guild_members_voice_count, guild_snowflake=channel.guild.id, highest_role=highest_role, is_modification=is_modification, target_member_snowflake=member.id, reason=reason)

    @classmethod
    def build_history_embeds(
        cls,
        alias,
        channel,
        duration,
        highest_role,
        is_channel_scope,
        is_modification,
        member,
        message,
        reason
    ):
        if duration != "permanent":
            if duration is None:
                duration_info = f'**Type:** {alias.alias_type.capitalize()}\n\n**Expires:** Never'
            else:
                duration_info = f'**Type:** {alias.alias_type.capitalize()}\n\n**Expires:** {duration}'
            if is_modification:
                color, duration_type = 0xFF6B35, '⏰ Modified'
            else:
                color, duration_type = 0xFF8C00, '⏱️ Temporary'
        else:
            color, duration_type = 0xDC143C, '♾️ Permanent'
            duration_info = f'**Type:** {alias.alias_type.capitalize()}\n**Expires:** {duration}'
        if alias.alias_type == 'ban':
            if is_modification:
                title = '🔄 Ban Modified'
            else:
                title = '🔨 User Banned'
            action = 'banned'
        elif alias.alias_type == 'flag':
            if is_modification:
                title = '🔄 Flag Modified'
            else:
                title = '🚩 User Flagged'
            action = 'flagged'
        elif alias.alias_type == 'text_mute':
            if is_modification:
                title = '🔄 Text Mute Modified'
            else:
                title = '📝 User Text Muted'
            action = 'text muted'
        elif alias.alias_type == 'voice_mute':
            if is_modification:
                title = '🔄 Voice Mute Modified'
            else:
                title = '🎙️ User Voice Muted'
            action = 'voice muted'
        elif alias.alias_type == 'unban':
            title = '⏮️ User Unbanned'
            action = 'unbanned'
        elif alias.alias_type == 'unflag':
            title = '⏮️ User Unflagged'
            action = 'unflagged'
        elif alias.alias_type == 'untext_mute':
            title = '⏮️ User Text-Mute Removed'
            action = 'untext-muted'
        elif alias.alias_type == 'unvoice_mute':
            title = '⏮️ User Voice-Mute Removed'
            action = 'unvoice-muted'
        else:
            title = None
            action = None
        embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
        if not action:
            embed_user.description = None
        else:
            embed_user.description = f"**Target:** {member.mention} {action} in {channel.mention}"
        embed_user.set_thumbnail(url=message.author.display_avatar.url)
        user_priority = f"**Display Name:** {member.display_name}\n**Username:** @{member.name}\n**User ID:** `{member.id}`\n**Account Age:** <t:{int(member.created_at.timestamp())}:R>\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
        embed_user.add_field(name='👤 Target User', value=user_priority, inline=False)
        exec_priority = f"**Executor:** {message.author.display_name} (@{message.author.name})\n**Executor ID:** `{message.author.id}`\n**Top Role:** {highest_role}"
        embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
        ctx_info = f"**Original Message ID:** `{message.id}`\n**Message Link:** [Jump to Message]({message.jump_url})\n**Command Channel:** {message.channel.mention}\n**Command Used:** `{alias.alias_name}`"
        embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
        embed_user.add_field(name=f'**Type:** {duration_type}', value=duration_info, inline=False)
        embed_user.add_field(name='📝 Reason', value=f'```{reason}```', inline=False)
        embed_user.set_footer(text=f"Ref: {member.id}-{channel.id} | Msg: {message.id}", icon_url=message.guild.icon.url)
        embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
        embed_duration.add_field(name=f'{duration_type}', value=duration_info, inline=False)
        action_details = f"**Was in Channel:** {is_channel_scope}\n**Modification:** {is_modification}\n**Server:** {message.guild.name}"
        embed_duration.add_field(name='⚙️ Action Details', value=action_details, inline=True)
        channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name}"
        embed_duration.add_field(name='📍 Channel Info', value=channel_basic, inline=True)
        embeds = [embed_user, embed_duration]
        if reason:
            reason_chunks = [reason[i:i+1000] for i in range(0, len(reason), 1000)]
            if len(reason_chunks) > 1:
                for i, chunk in enumerate(reason_chunks):
                    reason_embed = discord.Embed(title=f"{title} - Reason (cont.)", color=color, timestamp=datetime.now(timezone.utc))
                    reason_embed.add_field(name=f'📝 Reason (Part {i+1})', value=f'```{chunk}```', inline=False)
                    embeds.append(reason_embed)
        return embeds

    @classmethod
    async def delete_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                DELETE FROM history
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)

    @classmethod
    async def fetch_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, enabled, guild_snowflake, snowflakes, entry_type, updated_at
                FROM history
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
        history = []
        if rows:
            for row in rows:
                history.append(History(channel_snowflake=channel_snowflake, enabled=row['enabled'], guild_snowflake=guild_snowflake, snowflakes=row['snowflakes'], entry_type=row['entry_type']))
            return history

    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, enabled, guild_snowflake, snowflakes, entry_type, updated_at
                FROM history
                WHERE guild_snowflake=$1
            ''', guild_snowflake)
        history = []
        if rows:
            for row in rows:
                history.append(History(channel_snowflake=row['channel_snowflake'], enabled=row['enabled'], guild_snowflake=guild_snowflake, snowflakes=row['snowflakes'], entry_type=row['entry_type']))
            return history

    @classmethod
    async def update_by_channel_guild_and_type(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], snowflakes: list[int|None], entry_type: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE history
                SET snowflakes=$3, entry_type=$4, enabled=TRUE
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake, snowflakes, entry_type)

    @classmethod
    async def update_by_channel_enabled_and_guild(cls, channel_snowflake: Optional[int], enabled: Optional[bool], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE history
                SET enabled=$2
                WHERE channel_snowflake=$1 AND guild_snowflake=$3
            ''', channel_snowflake, enabled, guild_snowflake)

    async def create(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO history (channel_snowflake, enabled, guild_snowflake, snowflakes, entry_type)
                VALUES ($1, TRUE, $2, $3, $4)
            ''', self.channel_snowflake, self.guild_snowflake, self.snowflakes, self.entry_type)

    @classmethod
    async def fetch_all(cls):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, created_at, enabled, guild_snowflake, snowflakes, entry_type, updated_at
                FROM history
            ''')
        history = []
        if rows:
            for row in rows:
                history.append(
                    History(
                        channel_snowflake=row['channel_snowflake'],
                        enabled=row['enabled'],
                        guild_snowflake=row['guild_snowflake'],
                        snowflakes=row['snowflakes'],
                        entry_type=row['entry_type']
                    )
                )
        return history

    @classmethod
    async def save(
        cls,
        action_type: Optional[str],
        channel_members_voice_count: Optional[int],
        channel_snowflake: Optional[int],
        executor_member_snowflake: Optional[int],
        expires_at: Optional[datetime],
        guild_members_offline_and_online_member_count: Optional[int],
        guild_members_online_count: Optional[int],
        guild_members_voice_count: Optional[int],
        guild_snowflake: Optional[int],
        highest_role: Optional[str],
        is_modification: bool,
        target_member_snowflake: Optional[int],
        reason: Optional[str]
    ):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                    INSERT INTO moderation_logs (action_type, channel_members_voice_count, channel_snowflake, executor_member_snowflake, expires_at, guild_members_offline_and_online_member_count, guild_members_online_count, guild_members_voice_count, guild_snowflake, highest_role, is_modification, target_member_snowflake, reason)
                    VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13)
                ''', action_type, channel_members_voice_count, channel_snowflake, executor_member_snowflake, expires_at, guild_members_offline_and_online_member_count, guild_members_online_count, guild_members_voice_count, guild_snowflake, highest_role, is_modification, target_member_snowflake, reason)

    @classmethod
    async def save_entry(
        cls,
        ctx_interaction_or_message,
        action_type: Optional[str],
        channel_snowflake: Optional[int],
        duration: Optional[str],
        highest_role: Optional[str],
        is_modification: bool,
        member_snowflake: Optional[int],
        reason: Optional[str]
    ):
        bot = DiscordBot.get_instance()
        author_snowflake = None
        expires_at = None
        if isinstance(ctx_interaction_or_message, (commands.Context, discord.Message)):
            author_snowflake = ctx_interaction_or_message.author.id
        else:
            author_snowflake = ctx_interaction_or_message.user.id
        channel = ctx_interaction_or_message.channel
        if isinstance(duration, DurationObject):
            expires_at = datetime.now(timezone.utc) + duration.to_timedelta()
        elif duration is not None:
            expires_at = datetime.now(timezone.utc) + DurationObject(duration).to_timedelta()
        else:
            expires_at = datetime.now(timezone.utc) + DurationObject(0).to_timedelta()
        channel_members_voice_count = sum(len(channel.members) for channel in channel.guild.voice_channels)
        guild_members_offline_and_online_member_count = sum(1 for member in channel.guild.members if not member.bot)
        guild_members_online_count = sum(1 for member in channel.guild.members if not member.bot and member.status != discord.Status.offline)
        guild_members_voice_count = sum(len([member for member in channel.members if not member.bot]) for channel in channel.guild.voice_channels)
        await cls.save(action_type=action_type, channel_members_voice_count=channel_members_voice_count, channel_snowflake=channel_snowflake, executor_member_snowflake=author_snowflake, expires_at=expires_at, guild_members_offline_and_online_member_count=guild_members_offline_and_online_member_count, guild_members_online_count=guild_members_online_count, guild_members_voice_count=guild_members_voice_count, guild_snowflake=channel.guild.id, highest_role=highest_role, is_modification=is_modification, target_member_snowflake=member_snowflake, reason=reason)
