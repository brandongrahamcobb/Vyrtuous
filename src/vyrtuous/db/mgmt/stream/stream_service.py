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

import discord

from vyrtuous.base.service import Service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.commands.messaging.message_service import PaginatorService
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.mgmt.stream.stream import Stream
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.data import Data
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_channels,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji


class StreamService(Service):

    lines, pages = [], []

    @classmethod
    async def send_entry(
        cls,
        channel_snowflake: int,
        identifier: str,
        member: discord.Member,
        message: discord.Message | None = None,
        is_channel_scope: bool = False,
        is_modification: bool = False,
        reason: str = "No reason provided",
        duration: str | None = None,
    ):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(channel_snowflake)
        expires_at = None
        streaming = await Stream.select(singular=False)
        if message:
            executor_highest_role = await PermissionService.resolve_highest_role_at_all(
                member_snowflake=int(message.author.id),
            )
        else:
            executor_highest_role = "Role undetermined"
        target_highest_role = await PermissionService.resolve_highest_role_at_all(
            member_snowflake=int(member.id),
        )
        if message:
            for stream in streaming:
                channel_obj = bot.get_channel(stream.channel_snowflake)
                if channel_obj:
                    perms = channel_obj.permissions_for(channel_obj.guild.me)
                    if perms.send_messages and not channel_obj.guild.me.is_timed_out():
                        pages = cls.build_streaming_embeds(
                            channel=channel_obj,
                            duration=duration,
                            executor_highest_role=executor_highest_role,
                            identifier=identifier,
                            is_channel_scope=is_channel_scope,
                            is_modification=is_modification,
                            member=member,
                            message=message,
                            reason=reason,
                            target_highest_role=target_highest_role,
                        )
                        paginator = PaginatorService(bot, channel_obj, pages)
                        await paginator.start()
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
            channel_members_voice_count=channel_members_voice_count,
            channel_snowflake=int(channel.id),
            executor_member_snowflake=int(message.author.id) if message else None,
            expires_at=expires_at,
            guild_members_offline_and_online_member_count=guild_members_offline_and_online_member_count,
            guild_members_online_count=guild_members_online_count,
            guild_members_voice_count=guild_members_voice_count,
            guild_snowflake=int(channel.guild.id),
            executor_highest_role=executor_highest_role,
            identifier=identifier,
            is_modification=is_modification,
            target_member_snowflake=int(member.id),
            target_highest_role=target_highest_role,
            reason=reason,
        )

    @classmethod
    def build_streaming_embeds(
        cls,
        channel: discord.abc.GuildChannel,
        executor_highest_role: str,
        identifier: str,
        member: discord.Member,
        target_highest_role: str,
        duration: str | None = None,
        is_channel_scope: bool = False,
        is_modification: bool = False,
        message: discord.Message | None = None,
        reason: str = "No reason provided",
    ):
        if duration != "permanent":
            if duration is None:
                duration_info = (
                    f"**Type:** {identifier.capitalize()}\n**Expires:** Never"
                )
            else:
                duration_info = (
                    f"**Type:** {identifier.capitalize()}\n**Expires:** {duration}"
                )
            if is_modification:
                color, duration_type = 0xFF6B35, "⏰ Modified"
            else:
                color, duration_type = 0xFF8C00, "⏱️ Temporary"
        else:
            color, duration_type = 0xDC143C, "♾️ Permanent"
            duration_info = (
                f"**Type:** {identifier.capitalize()}\n**Expires:** {duration}"
            )
        if identifier in ("ban", "unban"):
            if is_modification:
                title = "🔄 User Unbanned"
            else:
                title = "🔨 User Banned"
            action = "banned"
        elif identifier in ("flag", "unflag"):
            if is_modification:
                title = "🔄 User Unflagged"
            else:
                title = "🚩 User Flaged"
            action = "flagged"
        elif identifier in ("tmute", "untmute"):
            if is_modification:
                title = "🔄 User Unmuted"
            else:
                title = "📝 User Text Muted"
            action = "text muted"
        elif identifier in ("vmute", "unvmute"):
            if is_modification:
                title = "🔄 User Unmuted"
            else:
                title = "🎙️ User Voice Muted"
            action = "voice muted"
        else:
            title = None
            action = None
        embed_user = discord.Embed(
            title=f"{title} - User Identity",
            color=color,
            timestamp=datetime.now(timezone.utc),
        )
        if not action:
            embed_user.description = None
        else:
            embed_user.description = (
                f"**Target:** {member.mention} {action} in {channel.mention}"
            )
        embed_user.set_thumbnail(url=message.author.display_avatar.url)
        user_priority = f"**Display Name:** {member.display_name}\n**Username:** @{member.name}\n**User ID:** `{member.id}`\n**Account Age:** <t:{int(member.created_at.timestamp())}:R>\n**Server Join:** <t:{int(member.joined_at.timestamp())}:R>\n**Top Role:** {target_highest_role}"
        embed_user.add_field(name="👤 Target User", value=user_priority, inline=False)
        exec_priority = f"**Executor:** {message.author.display_name} (@{message.author.name})\n**Executor ID:** `{message.author.id}`\n**Top Role:** {executor_highest_role}"
        embed_user.add_field(name="👮‍♂️ Executed By", value=exec_priority, inline=True)
        ctx_info = f"**Original Message ID:** `{message.id}`\n**Message Link:** [Jump to Message]({message.jump_url})\n**Command Channel:** {message.channel.mention}\n**Type:** `{identifier}`"
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
        details = f"**Was in Channel:** {is_channel_scope}\n**Modification:** {is_modification}\n**Server:** {message.guild.name}"
        embed_duration.add_field(name="⚙️ Action Details", value=details, inline=True)
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
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        streaming = await Stream.select(singular=False, **where_kwargs)
        for stream in streaming:
            dictionary.setdefault(stream.guild_snowflake, {"channels": {}})
            dictionary[stream.guild_snowflake]["channels"][stream.channel_snowflake] = {
                "enabled": stream.enabled,
                "entry_type": stream.entry_type,
                "snowflakes": stream.snowflakes,
            }
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_channels:
                StreamService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
            if skipped_guilds:
                StreamService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Streaming Routes"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await StreamService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            stream_n = 0
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, entry in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                status = "\u2705" if entry["enabled"] else "\u26d4"
                StreamService.lines.append(
                    f"{status}**Channel:** {channel.mention}\n**Type:** {entry['entry_type']}"
                )
                if isinstance(
                    object_dict.get("object", None), discord.abc.GuildChannel
                ):
                    StreamService.lines.append(f"**Snowflakes:** {entry['snowflakes']}")
                field_count += 1
                stream_n += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(StreamService.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, StreamService.pages, title, guild.name)
                    StreamService.lines = []
                    field_count = 0
            if StreamService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(StreamService.lines),
                    inline=False,
                )
            StreamService.pages.append(embed)
            StreamService.pages[0].description = f'{guild.name} **({stream_n})**'
        return StreamService.pages

    @classmethod
    async def modify_stream(
        cls,
        action,
        channel_dict,
        channel_mentions,
        default_kwargs,
        entry_type,
        failed_snowflakes,
        resolved_channels,
    ):
        updated_kwargs = default_kwargs.copy()
        guild_snowflake = updated_kwargs.get("guild_snowflake", None)
        channel_kwargs = channel_dict.get("columns", None)
        where_kwargs = {
            "channel_snowflake": channel_kwargs["channel_snowflake"],
            "entry_type": entry_type,
            "guild_snowflake": int(guild_snowflake),
        }
        if action is None and entry_type is None:
            stream = await Stream.select(singular=True, **where_kwargs)
            enabled = not stream.enabled
            action = "enabled" if enabled else "disabled"
            set_kwargs = {"enabled": enabled}
            await Stream.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        if action and entry_type:
            match action.lower():
                case "create" | "modify":
                    resolved_channels = []
                    if action.lower() == "create":
                        stream = Stream(
                            **where_kwargs,
                            enabled=True,
                            snowflakes=resolved_channels,
                        )
                        await stream.create()
                        action = "created"
                    else:
                        set_kwargs = {"snowflakes": resolved_channels}
                        await Stream.update(
                            set_kwargs=set_kwargs, where_kwargs=where_kwargs
                        )
                        action = "modified"
                case "delete":
                    await Stream.delete(**where_kwargs)
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
