
''' statistics.py A utility module for managing and sending statistical messages in the Vyrtuous Discord bot.
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
from vyrtuous.service.discord_message_service import Paginator

import discord

class Statistics:
        
    ACTION_TYPES = ['create', 'delete', 'modify']
    STATISTIC_TYPES = ['channel', 'general', 'member']

    def __init__(self, channel_snowflake: Optional[int], enabled: Optional[bool], guild_snowflake: Optional[int], snowflakes: list[int|None], statistic_type: Optional[str]):
        self.action: Optional[str]
        self.channel_snowflake = channel_snowflake
        self.enabled = enabled
        self.guild_snowflake = guild_snowflake
        self.snowflakes = snowflakes
        self.statistic_type = statistic_type

    @classmethod
    async def send_statistic(
        cls,
        message,
        moderation_type: Optional[str],
        member: discord.Member,
        channel: Optional[discord.VoiceChannel],
        duration_display: Optional[str],
        reason: Optional[str],
        expires_at: Optional[datetime],
        command_used: Optional[str],
        was_in_channel: bool = False,
        is_modification: bool = False,
        highest_role: Optional[str] = 'Everyone'
    ):
        bot = DiscordBot.get_instance()
        statistics = await Statistics.fetch_all()
        for statistic in statistics:
            if message.guild.id == statistic.channel_snowflake:
                channel = bot.get_channel(statistic.channel_snowflake)
                statistic_type = statistic.statistic_type
                snowflakes = statistic.snowflakes
                pages = cls.build_statistic_embeds(
                        message=message, moderation_type=moderation_type,member=member, channel=channel, duration_display=duration_display,
                        reason=reason, executor=message.author, expires_at=expires_at,
                        command_used=command_used, was_in_channel=was_in_channel,
                        is_modification=is_modification, highest_role=highest_role
                    )
                paginator = Paginator(bot, channel, pages)
                await paginator.start()

    @classmethod
    def build_statistic_embeds(cls, message: discord.Message, moderation_type: Optional[str], member: discord.Member, channel: discord.VoiceChannel, duration_display: Optional[str], reason: Optional[str], executor: discord.Member, expires_at: Optional[datetime], command_used: Optional[str], was_in_channel: bool = False, is_modification: bool = False, highest_role: Optional[str] = ''):
        guild = message.guild
        if expires_at is not None:
            time_left = expires_at - datetime.now(timezone.utc)
            hours_left = round(time_left.total_seconds() / 3600, 1)
            days_left = time_left.days
            duration_info = f'**Type:** {moderation_type}\n**Duration:**\n{duration_display}\n**Expires:** <t:{int(expires_at.timestamp())}:F>\n**Time Left:** '
            duration_info += f'{days_left}d {hours_left % 24:.1f}h' if days_left > 0 else f'{hours_left}h'
            if is_modification:
                color, duration_type = 0xFF6B35, '⏰ Modified'
            else:
                color, duration_type = 0xFF8C00, '⏱️ Temporary'
        else:
            color, duration_type = 0xDC143C, '♾️ Permanent'
            duration_info = f'**Type:** {moderation_type}\n**Duration:** {duration_display}\n**Status:** Permanent'
            
        if moderation_type and moderation_type.lower() == 'ban':
            if is_modification:
                title = '🔄 Ban Modified'
            else:
                title = '🔨 User Banned'
            action = 'banned'
        elif moderation_type and moderation_type.lower() == 'voice_mute':
            if is_modification:
                title = '🔄 Voice Mute Modified'
            else:
                title = '🎙️ User Voice Muted'
            action = 'voice muted'
        elif moderation_type and moderation_type.lower() == 'text_mute':
            if is_modification:
                title = '🔄 Text Mute Modified'
            else:
                title = '📝 User Text Muted'
            action = 'text muted'
        else:
            title = None
            action = None
            
        embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
        if not action:
            embed_user.description = None
        else:
            embed_user.description = f"**Target:** {member.mention} {action} in {channel.mention}"
        embed_user.set_thumbnail(url=message.author.display_avatar.url)
        user_priority = f"**Display Name:** {member.display_name}\n**Username:** @{member.name}\n**User ID:** `{member.id}`\n**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
        if member.joined_at:
            user_priority += f"\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
        embed_user.add_field(name='👤 Target User', value=user_priority, inline=False)
        exec_priority = f"**Moderator:** {message.author.display_name} (@{message.author.name})\n**Mod ID:** `{message.author.id}`\n**Top Role:** {highest_role or message.author.top_role.mention}"
        embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
        ctx_info = f"**Original Message ID:** `{message.id}`\n**Message Link:** [Jump to Message]({message.jump_url})\n**Command Channel:** {message.channel.mention}\n**Command Used:** `{command_used}`"
        embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
        embed_user.add_field(name=f'**Type:** {duration_type}', value=duration_info, inline=False)
        embed_user.add_field(name='📝 Reason', value=f'```{reason if reason else "No reason provided"}```', inline=False)
        embed_user.set_footer(text=f"Ref: {member.id}-{channel.id} | Msg: {message.id}", icon_url=guild.icon.url if guild and guild.icon else None)
        embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
        embed_duration.add_field(name=f'{duration_type}', value=duration_info, inline=False)
        action_details = f"**Was in Channel:** {'✅ Yes' if was_in_channel else '\U0001F6AB No'}\n**Action Type:** {'Modification' if is_modification else 'New'}\n**Server:** {guild.name} (`{guild.id}`)"
        embed_duration.add_field(name='⚙️ Action Details', value=action_details, inline=True)
        channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name if channel.category else 'None'}"
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
                DELETE FROM statistic_channels
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)

    @classmethod
    async def fetch_by_channel_and_guild(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, enabled, guild_snowflake, snowflakes, statistic_type
                FROM statistic_channels
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake)
        statistics = []
        if rows:
            for row in rows:
                statistics.append(Statistics(channel_snowflake=channel_snowflake, enabled=row['enabled'], guild_snowflake=guild_snowflake, snowflakes=row['snowflakes'], statistic_type=row['statistic_type']))
            return statistics

    @classmethod
    async def fetch_by_guild(cls, guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch('''
                SELECT channel_snowflake, enabled, guild_snowflake, snowflakes, statistic_type
                FROM statistic_channels
                WHERE guild_snowflake=$1
            ''', guild_snowflake)
        statistics = []
        if rows:
            for row in rows:
                statistics.append(Statistics(channel_snowflake=row['channel_snowflake'], enabled=row['enabled'], guild_snowflake=guild_snowflake, snowflakes=row['snowflakes'], statistic_type=row['statistic_type']))
            return statistics

    @classmethod
    async def update_by_channel_guild_and_type(cls, channel_snowflake: Optional[int], guild_snowflake: Optional[int], snowflakes: list[int|None], statistic_type: Optional[str]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE statistic_channels
                SET snowflakes=$3, statistic_type=$4, enabled=TRUE
                WHERE channel_snowflake=$1 AND guild_snowflake=$2
            ''', channel_snowflake, guild_snowflake, snowflakes, statistic_type)

    @classmethod
    async def update_by_channel_enabled_and_guild(cls, channel_snowflake: Optional[int], enabled: Optional[bool], guild_snowflake: Optional[int]):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                UPDATE statistic_channels
                SET enabled=$2
                WHERE channel_snowflake=$1 AND guild_snowflake=$3
            ''', channel_snowflake, enabled, guild_snowflake)

    async def grant(self):
        bot = DiscordBot.get_instance()
        async with bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO statistic_channels (channel_snowflake, enabled, guild_snowflake, snowflakes, statistic_type)
                VALUES ($1, TRUE, $2, $3, $4)
            ''', self.channel_snowflake, self.guild_snowflake, self.snowflakes, self.statistic_type)

    @property
    def action(self):
        self._action

    @action.setter
    def action(self, action: Optional[str]):
        if action not in self.ACTION_TYPES:
            raise ValueError('Invalid action.')
        self._action = action

    @property
    def statistic_type(self):
        return self._statistic_type

    @statistic_type.setter
    def statistic_type(self, statistic_type: Optional[str]):
        if statistic_type not in self.STATISTIC_TYPES:
            raise ValueError('Invalid statistic type.')
        self._statistic_type = statistic_type