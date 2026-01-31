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
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.data import Data
from vyrtuous.utils.highest_role import resolve_highest_role
from vyrtuous.fields.duration import DurationObject
from vyrtuous.service.message_service import PaginatorService
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_guild_dictionary,
    flush_page,
)


class Streaming(DatabaseFactory):

    ACTION_TYPES = ["create", "delete", "modify"]
    ENTRY_TYPES = ["all", "channel", "member"]
    ACT = None
    CATEGORY = "stream"
    PLURAL = "Streaming Channels"
    SCOPES = ["channels"]
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
        if action not in Streaming.ACTION_TYPES:
            raise ValueError("Invalid action.")
        self._action = action

    @property
    def entry_type(self):
        return self._entry_type

    @entry_type.setter
    def entry_type(self, entry_type: str):
        if entry_type not in Streaming.ENTRY_TYPES:
            raise ValueError("Invalid entry type.")
        self._entry_type = entry_type

    @classmethod
    async def send_entry(
        cls,
        alias,
        channel_snowflake: int,
        duration: str,
        is_channel_scope: bool,
        is_modification: bool,
        member: discord.Member,
        message,
        reason: str,
    ):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(channel_snowflake)
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
                    paginator = PaginatorService(bot, channel_obj, pages)
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
            action_type=alias.category,
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
                    f"**Type:** {alias.category.capitalize()}\n**Expires:** Never"
                )
            else:
                duration_info = (
                    f"**Type:** {alias.category.capitalize()}\n**Expires:** {duration}"
                )
            if is_modification:
                color, duration_type = 0xFF6B35, "⏰ Modified"
            else:
                color, duration_type = 0xFF8C00, "⏱️ Temporary"
        else:
            color, duration_type = 0xDC143C, "♾️ Permanent"
            duration_info = (
                f"**Type:** {alias.category.capitalize()}\n**Expires:** {duration}"
            )
        if alias.category == "ban":
            if is_modification:
                title = "🔄 Ban Modified"
            else:
                title = "🔨 User Banned"
            action = "banned"
        elif alias.category == "flag":
            if is_modification:
                title = "🔄 Flag Modified"
            else:
                title = "🚩 User Flagged"
            action = "flagged"
        elif alias.category == "tmute":
            if is_modification:
                title = "🔄 Text Mute Modified"
            else:
                title = "📝 User Text Muted"
            action = "text muted"
        elif alias.category == "vmute":
            if is_modification:
                title = "🔄 Voice Mute Modified"
            else:
                title = "🎙️ User Voice Muted"
            action = "voice muted"
        elif alias.category == "unban":
            title = "⏮️ User Unbanned"
            action = "unbanned"
        elif alias.category == "unflag":
            title = "⏮️ User Unflagged"
            action = "unflagged"
        elif alias.category == "untmute":
            title = "⏮️ User Text-Mute Removed"
            action = "untext-muted"
        elif alias.category == "unvmute":
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

    @classmethod
    async def build_dictionary(cls, where_kwargs):
        dictionary = {}
        streaming = await Streaming.select(**where_kwargs)
        for stream in streaming:
            dictionary.setdefault(stream.guild_snowflake, {"channels": {}})
            dictionary[stream.guild_snowflake]["channels"][
                stream.channel_snowflake
            ] = {
                "enabled": stream.enabled,
                "entry_type": stream.entry_type,
                "snowflakes": stream.snowflakes,
            }
        return dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        title = f"{get_random_emoji()} {Streaming.PLURAL}"

        where_kwargs = object_dict.get("columns", None)

        guild_dictionary = await Streaming.build_dictionary(where_kwargs=where_kwargs)

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, entry in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                status = "\u2705" if entry["enabled"] else "\u26d4"
                lines.append(
                    f"{status}**Channel:** {channel.mention}\n**Type:** {entry['entry_type']}"
                )
                if isinstance(
                    object_dict.get("object", None), discord.abc.GuildChannel
                ):
                    lines.append(f"**Snowflakes:** {entry['snowflakes']}")
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
                )
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
        return pages

    @classmethod
    async def modify_stream(
        cls,
        action,
        channel_dict,
        channel_mentions,
        entry_type,
        failed_snowflakes,
        resolved_channels,
        snowflake_kwargs,
    ):
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        channel_kwargs = channel_dict.get("columns", None)
        where_kwargs = {
            "channel_snowflake": channel_kwargs["channel_snowflake"],
            "entry_type": entry_type,
            "guild_snowflake": guild_snowflake,
        }
        if action is None and entry_type is None:
            stream = await Streaming.select(**channel_kwargs, singular=True)
            enabled = not stream[0].enabled
            action = "enabled" if enabled else "disabled"
            set_kwargs = {"enabled": enabled}
            await Streaming.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        if action and entry_type:
            match action.lower():
                case "create" | "modify":
                    resolved_channels = []
                    if action.lower() == "create":
                        stream = Streaming(
                            **channel_kwargs,
                            enabled=True,
                            entry_type=entry_type,
                            snowflakes=resolved_channels,
                        )
                        await stream.create()
                        action = "created"
                    else:
                        set_kwargs = {"snowflakes": resolved_channels}
                        await Streaming.update(
                            set_kwargs=set_kwargs, where_kwargs=where_kwargs
                        )
                        action = "modified"
                case "delete":
                    await Streaming.delete(**channel_kwargs)
                    action = "deleted"
        embed = discord.Embed(
            title=f"{get_random_emoji()} Tracking {action.capitalize()} for {channel_dict.get("mention", None)}",
            color=0x00FF00,
        )
        if channel_mentions:
            embed.add_field(
                name="Processed Channels",
                value=", ".join(channel_mentions),
                inline=False,
            )
        if failed_snowflakes:
            embed.add_field(
                name="Failed IDs",
                value=", ".join(str(s) for s in failed_snowflakes),
                inline=False,
            )
        return [embed]
