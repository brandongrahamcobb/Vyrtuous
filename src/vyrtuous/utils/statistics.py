
from datetime import datetime
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.discord_message_service import ChannelPaginator, UserPaginator

import discord

class Statistics:

    def __init__(self):
        self.bot = DiscordBot.get_instance()
        self.statistic_channels: dict[int, list[dict]] = {}

    async def load_channels(self):
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch('SELECT guild_id, channel_id, type, snowflakes, enabled FROM statistic_channels;')
            statistic_channels: dict[int, list[dict]] = {}
            for r in rows:
                statistic_channels.setdefault(r['guild_id'], []).append({
                    "channel_id": r['channel_id'],
                    "type": r['type'] or "general",
                    "snowflakes": r['snowflakes'] or [],
                    "enabled": r['enabled']
                })
            self.statistic_channels = statistic_channels

    async def send_statistic(
        self,
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
        guild_id = message.guild.id
        guild = message.guild
        if guild_id not in self.statistic_channels:
            print(f"No statistic channels for guild {guild_id}")
            return
        for entry in self.statistic_channels[guild_id]:
            log_channel = self.bot.get_channel(entry["channel_id"])
            if not log_channel:
                continue
            log_type = entry.get("type")
            snowflakes = entry.get("snowflakes") or []
            if log_type == "member" and member.id not in snowflakes:
                continue
            if log_type == "channel" and channel and channel.id not in snowflakes:
                continue
            pages = self.build_statistic_embeds(
                    message=message, member=member, channel=channel, duration_display=duration_display,
                    reason=reason, executor=executor, expires_at=expires_at,
                    command_used=command_used, was_in_channel=was_in_channel,
                    is_modification=is_modification, highest_role=highest_role
                )
            paginator = ChannelPaginator(self.bot, log_channel, pages)
            await paginator.start()

    def build_statistic_embeds(self, message: discord.Message, member: discord.Member, channel: discord.VoiceChannel, duration_display: Optional[str], reason: Optional[str], executor: discord.Member, expires_at: Optional[datetime], command_used: Optional[str], was_in_channel: bool = False, is_modification: bool = False, highest_role: Optional[str] = '', moderation_type: Optional[str] = None):
        guild = message.guild
        time_left = expires_at - datetime.now(timezone.utc)
        hours_left = round(time_left.total_seconds() / 3600, 1)
        days_left = time_left.days
        duration_info = f'**Type:** {moderation_type}\n**Duration:** {duration_display}\n**Expires:** <t:{int(expires_at.timestamp())}:F>\n**Time Left:** '
        duration_info += f'{days_left}d {hours_left % 24:.1f}h' if days_left > 0 else f'{hours_left}h'
        if expires_at is None:
            color, duration_type, duration_emoji = 0xDC143C, '🔒 Permanent', '♾️'
            duration_info = f'**Type:** {moderation_type}\n**Duration:** {duration_display}\n**Status:** Permanent'
        elif is_modification:
            color, duration_type, duration_emoji = 0xFF6B35, '⏰ Modified', '📅'
        else:
            color, duration_type, duration_emoji = 0xFF8C00, '⏱️ Temporary', '⏰'
            
        if moderation_type.lower() == 'ban':
            if is_modification:
                title = '🔄 Ban Modified'
            else:
                title = '🔨 User Banned'
            action = 'banned'
        elif moderation_type.lower() == 'voice_mute':
            if is_modification:
                title = '🔄 Voice Mute Modified'
            else:
                title = '🎙️ User Voice Muted'
            action = 'voice muted'
        elif moderation_type.lower() == 'text_mute':
            if is_modification:
                title = '🔄 Text Mute Modified'
            else:
                title = '🔇 User Text Muted'
            action = 'text muted'
        else:
            title = None
            embed_user.description = None
            action = None
            
        embed_user = discord.Embed(title=f"{title} - User Identity", color=color, timestamp=datetime.now(timezone.utc))
        embed_user.description = f"**Target:** {member.mention} {action} in {channel.mention}"
        embed_user.set_thumbnail(url=message.author.display_avatar.url)
        user_priority = f"**Display Name:** {member.display_name}\n**Username:** @{member.name}\n**User ID:** `{member.id}`\n**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
        if member.joined_at:
            user_priority += f"\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
        embed_user.add_field(name='👤 Target User (Priority Info)', value=user_priority, inline=False)
        exec_priority = f"**Moderator:** {message.author.display_name} (@{message.author.name})\n**Mod ID:** `{message.author.id}`\n**Top Role:** {highest_role or message.author.top_role.mention}"
        embed_user.add_field(name='👮‍♂️ Executed By', value=exec_priority, inline=True)
        ctx_info = f"**Original Message ID:** `{message.id}`\n**Message Link:** [Jump to Message]({message.jump_url})\n**Command Channel:** {message.channel.mention}\n**Command Used:** `{command_used}`"
        embed_user.add_field(name='📱 Command Context', value=ctx_info, inline=True)
        embed_user.add_field(name=f'{duration_emoji} Duration', value=duration_info, inline=False)
        embed_user.add_field(name='📝 Reason', value=f'```{reason if reason else "No reason provided"}```', inline=False)
        embed_user.set_footer(text=f"Ban Ref: {member.id}-{channel.id} | Msg: {message.id}", icon_url=guild.icon.url if guild and guild.icon else None)
        embed_duration = discord.Embed(title=f"{title} - Duration Info", color=color, timestamp=datetime.now(timezone.utc))
        embed_duration.add_field(name=f'{duration_emoji} Duration', value=duration_info, inline=False)
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
