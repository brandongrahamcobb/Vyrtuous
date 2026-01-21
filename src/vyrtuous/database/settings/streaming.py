"""streaming.py A utility module for managing and streaming of messages in the Vyrtuous Discord bot.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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
"""

from datetime import datetime, timezone
from typing import Optional

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.database_factory import DatabaseFactory
from vyrtuous.database.logs.data import Data
from vyrtuous.database.roles.role import resolve_highest_role
from vyrtuous.properties.duration import DurationObject
from vyrtuous.service.messaging.paginator_service import Paginator


class Streaming(DatabaseFactory):

    ACTION_TYPES = ["create", "delete", "modify"]
    ENTRY_TYPES = ["all", "channel", "member"]
    ACT = None
    PLURAL = "Streaming Channels"
    SINGULAR = "Streaming Channel"
    UNDO = None
    REQUIRED_INSTANTIATION_ARGS = [
        "channel_snowflake",
        "enabled",
        "entry_type",
        "guild_snowflake",
        "message_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "snowflakes",
        "updated_at",
    ]
    TABLE_NAME = "streaming"

    def __init__(
        self,
        channel_snowflake: int,
        enabled: Optional[bool],
        entry_type: str,
        guild_snowflake: int,
        created_at: Optional[datetime] = None,
        snowflakes: list[int] = None,
        updated_at: Optional[datetime] = None,
    ):
        self.action: str
        self.channel_snowflake = channel_snowflake
        self.created_at = created_at
        self.enabled = enabled
        self.entry_type = entry_type
        self.guild_snowflake = guild_snowflake
        self.snowflakes = snowflakes
        self.updated_at = updated_at

    @property
    def action(self):
        self._action

    @action.setter
    def action(self, action: str):
        if action not in self.ACTION_TYPES:
            raise ValueError("Invalid action.")
        self._action = action

    @property
    def entry_type(self):
        return self._entry_type

    @entry_type.setter
    def entry_type(self, entry_type: str):
        if entry_type not in self.ENTRY_TYPES:
            raise ValueError("Invalid entry type.")
        self._entry_type = entry_type

    @classmethod
    async def send_entry(
        cls,
        alias,
        channel: Optional[discord.VoiceChannel],
        duration: str,
        is_channel_scope: bool,
        is_modification: bool,
        member: discord.Member,
        message,
        reason: str,
    ):
        bot = DiscordBot.get_instance()
        author_snowflake = None
        expires_at = None
        streaming = await Streaming.select()
        highest_role = await resolve_highest_role(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
        )
        if message:
            for stream in streaming:
                channel_obj = bot.get_channel(stream.channel_snowflake)
                if channel_obj:
                    pages = cls.build_streaming_embeds(
                        alias=alias,
                        channel=channel_obj,
                        duration=duration,
                        highest_role=highest_role,
                        is_channel_scope=is_channel_scope,
                        is_modification=is_modification,
                        member=member,
                        message=message,
                        reason=reason,
                    )
                    paginator = Paginator(bot, channel_obj, pages)
                    await paginator.start()
            author_snowflake = message.author.id
        if isinstance(duration, DurationObject):
            expires_at = datetime.now(timezone.utc) + duration.to_timedelta()
        elif duration is not None:
            expires_at = (
                datetime.now(timezone.utc) + DurationObject(duration).to_timedelta()
            )
        else:
            expires_at = datetime.now(timezone.utc) + DurationObject(0).to_timedelta()
        channel_members_voice_count = len(channel.members)
        guild_members_offline_and_online_member_count = sum(
            1 for member in channel.guild.members if not member.bot
        )
        guild_members_online_count = sum(
            1
            for member in channel.guild.members
            if not member.bot and member.status != discord.Status.offline
        )
        guild_members_voice_count = sum(
            len([member for member in channel.members if not member.bot])
            for channel in channel.guild.voice_channels
        )
        await Data.save(
            action_type=alias.alias_type,
            channel_members_voice_count=channel_members_voice_count,
            channel_snowflake=channel.id,
            executor_member_snowflake=author_snowflake,
            expires_at=expires_at,
            guild_members_offline_and_online_member_count=guild_members_offline_and_online_member_count,
            guild_members_online_count=guild_members_online_count,
            guild_members_voice_count=guild_members_voice_count,
            guild_snowflake=channel.guild.id,
            highest_role=highest_role,
            is_modification=is_modification,
            target_member_snowflake=member.id,
            reason=reason,
        )

    @classmethod
    def build_streaming_embeds(
        cls,
        alias,
        channel,
        duration,
        highest_role,
        is_channel_scope,
        is_modification,
        member,
        message,
        reason,
    ):
        if duration != "permanent":
            if duration is None:
                duration_info = (
                    f"**Type:** {alias.alias_type.capitalize()}\n\n**Expires:** Never"
                )
            else:
                duration_info = f"**Type:** {alias.alias_type.capitalize()}\n\n**Expires:** {duration}"
            if is_modification:
                color, duration_type = 0xFF6B35, "⏰ Modified"
            else:
                color, duration_type = 0xFF8C00, "⏱️ Temporary"
        else:
            color, duration_type = 0xDC143C, "♾️ Permanent"
            duration_info = (
                f"**Type:** {alias.alias_type.capitalize()}\n**Expires:** {duration}"
            )
        if alias.alias_type == "ban":
            if is_modification:
                title = "🔄 Ban Modified"
            else:
                title = "🔨 User Banned"
            action = "banned"
        elif alias.alias_type == "flag":
            if is_modification:
                title = "🔄 Flag Modified"
            else:
                title = "🚩 User Flagged"
            action = "flagged"
        elif alias.alias_type == "tmute":
            if is_modification:
                title = "🔄 Text Mute Modified"
            else:
                title = "📝 User Text Muted"
            action = "text muted"
        elif alias.alias_type == "vmute":
            if is_modification:
                title = "🔄 Voice Mute Modified"
            else:
                title = "🎙️ User Voice Muted"
            action = "voice muted"
        elif alias.alias_type == "unban":
            title = "⏮️ User Unbanned"
            action = "unbanned"
        elif alias.alias_type == "unflag":
            title = "⏮️ User Unflagged"
            action = "unflagged"
        elif alias.alias_type == "untmute":
            title = "⏮️ User Text-Mute Removed"
            action = "untext-muted"
        elif alias.alias_type == "unvmute":
            title = "⏮️ User Voice-Mute Removed"
            action = "unvoice-muted"
        else:
            title = None
            action = None
        embed_user = discord.Embed(
            title=f"{title} - User Identity",
            color=color,
            timestamp=datetime.now(timezone.utc),
        )
        bot = DiscordBot.get_instance()
        channel_obj = bot.get_channel(alias.channel_snowflake)
        if not action:
            embed_user.description = None
        else:
            embed_user.description = (
                f"**Target:** {member.mention} {action} in {channel_obj.mention}"
            )
        embed_user.set_thumbnail(url=message.author.display_avatar.url)
        user_priority = f"**Display Name:** {member.display_name}\n**Username:** @{member.name}\n**User ID:** `{member.id}`\n**Account Age:** <t:{int(member.created_at.timestamp())}:R>\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
        embed_user.add_field(name="👤 Target User", value=user_priority, inline=False)
        exec_priority = f"**Executor:** {message.author.display_name} (@{message.author.name})\n**Executor ID:** `{message.author.id}`\n**Top Role:** {highest_role}"
        embed_user.add_field(name="👮‍♂️ Executed By", value=exec_priority, inline=True)
        ctx_info = f"**Original Message ID:** `{message.id}`\n**Message Link:** [Jump to Message]({message.jump_url})\n**Command Channel:** {message.channel.mention}\n**Command Used:** `{alias.alias_name}`"
        embed_user.add_field(name="📱 Command Context", value=ctx_info, inline=True)
        embed_user.add_field(
            name=f"**Type:** {duration_type}", value=duration_info, inline=False
        )
        embed_user.add_field(name="📝 Reason", value=f"```{reason}```", inline=False)
        embed_user.set_footer(
            text=f"Ref: {member.id}-{channel.id} | Msg: {message.id}",
            icon_url=message.guild.icon.url,
        )
        embed_duration = discord.Embed(
            title=f"{title} - Duration Info",
            color=color,
            timestamp=datetime.now(timezone.utc),
        )
        embed_duration.add_field(
            name=f"{duration_type}", value=duration_info, inline=False
        )
        action_details = f"**Was in Channel:** {is_channel_scope}\n**Modification:** {is_modification}\n**Server:** {message.guild.name}"
        embed_duration.add_field(
            name="⚙️ Action Details", value=action_details, inline=True
        )
        channel_basic = f"**Channel:** {channel.mention} (`{channel.id}`)\n**Category:** {channel.category.name}"
        embed_duration.add_field(
            name="📍 Channel Info", value=channel_basic, inline=True
        )
        embeds = [embed_user, embed_duration]
        if reason:
            reason_chunks = [reason[i : i + 1000] for i in range(0, len(reason), 1000)]
            if len(reason_chunks) > 1:
                for i, chunk in enumerate(reason_chunks):
                    reason_embed = discord.Embed(
                        title=f"{title} - Reason (cont.)",
                        color=color,
                        timestamp=datetime.now(timezone.utc),
                    )
                    reason_embed.add_field(
                        name=f"📝 Reason (Part {i+1})",
                        value=f"```{chunk}```",
                        inline=False,
                    )
                    embeds.append(reason_embed)
        return embeds
