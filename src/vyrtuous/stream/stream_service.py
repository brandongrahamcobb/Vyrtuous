"""!/bin/python3
streaming_service.py The purpose of this program is to extend Service service the stream command class.

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

from copy import copy
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Dict, List, Self, Union

import discord
from discord.ext import commands

from vyrtuous.stream.stream import Stream


@dataclass
class StreamDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, List[int]]]]] = field(
        default_factory=dict
    )
    skipped_channels: List[discord.Embed] = field(default_factory=list)
    skipped_guilds: List[discord.Embed] = field(default_factory=list)


class StreamEmbed(discord.Embed):
    def __init__(
        self,
        *,
        color=None,
        description=None,
        duration_builder=None,
        title=None,
        url=None,
    ):
        self.__action = None
        self.__duration_builder = duration_builder
        self.__timestamp = datetime.now(timezone.utc)
        super().__init__(
            color=color,
            description=description,
            timestamp=self.__timestamp,
            title=title,
            url=url,
        )

    def set_title(self, *, identifier, is_modification):
        if is_modification:
            self.color = 0xFF6B35
        else:
            self.color = 0xDC143C
        if identifier in ("ban", "unban"):
            if is_modification:
                self.title = "🔄 User Unbanned"
            else:
                self.title = "🔨 User Banned"
        elif identifier in ("flag", "unflag"):
            if is_modification:
                self.title = "🔄 User Unflagged"
            else:
                self.title = "🚩 User Flaged"
        elif identifier in ("tmute", "untmute"):
            if is_modification:
                self.title = "🔄 User Unmuted"
            else:
                self.title = "📝 User Text Muted"
        elif identifier in ("vmute", "unvmute"):
            if is_modification:
                self.title = "🔄 User Unmuted"
            else:
                self.title = "🎙️ User Voice Muted"
        return self

    def set_tn(self, *, url):
        self.set_thumbnail(url=url)
        return self

    def set_executor(self, *, author, highest_role) -> Self:
        fields = []
        if author:
            fields.append(f"**Executor:** {author.display_name} (@{author.name})")
            fields.append(f"**Executor id:** `{author.id}`")
        else:
            fields.append(f"**Executor:** Unknown")
            fields.append(f"**Executor id:** Unknown")
        fields.append(f"**Top Role:** {highest_role}")
        field = "\n".join(fields)
        self.add_field(name="👮‍♂️ Executed By", value=field, inline=True)
        return self

    def set_target(self, *, target, highest_role) -> Self:
        fields = []
        full_target = hasattr(target, "id")
        if full_target:
            fields.append(f"**Display Name:** {target.display_name}")
            fields.append(f"**Username:** @{target.name}")
            fields.append(f"**User ID:** `{target.id}`")
            fields.append(f"**Account Age:** <t:{int(target.created_at.timestamp())}:R>")
            fields.append(f"**Server Join:** <t:{int(target.joined_at.timestamp())}:R>")
        else:
            fields.append(f"**Display Name:** {target.get('name', None)}")
            fields.append(f"**User ID:** `{target.get('id', None)}`")
        fields.append(f"**Top Role:** {highest_role}")
        field = "\n".join(fields)
        self.add_field(name="👤 Target User", value=field, inline=False)
        return self

    def set_message_ctx(self, *, identifier, message) -> Self:
        fields = []
        fields.append(f"**Original Message ID:** `{message.id}`")
        fields.append(f"**Message Link:** [Jump to Message]({message.jump_url})")
        fields.append(f"**Command Channel:** {message.channel.mention}")
        fields.append(f"**Type:** `{identifier}`")
        field = "\n".join(fields)
        self.add_field(name="📱 Message Context", value=field, inline=True)
        return self

    def set_action(self, *, duration_value) -> Self:
        expiration = f"{self.__duration_builder.parse(duration_value).to_unix_ts()}"
        if duration_value is not None:
            dt = "⏱️ Temporary"
        else:
            dt = "♾️ Permanent"
            expiration = "Never"
        fields = []
        fields.append(f"**Expires:** {expiration}")
        self.add_field(name=f"**Type:** {dt}", value="\n".join(fields), inline=False)
        return self

    def set_reference(self, *, channel, target, message, source) -> Self:
        if hasattr(target, "id"):
            text = f"Ref: {target.id}-{channel.id} | Msg: {message.id if message else 'Hidden'}"
        else:
            text = f"Ref: {target.get("id", None)}-{channel.id} | Msg: {message.id if message else 'Hidden'}"
        icon_url = source.guild.icon.url
        self.set_footer(text=text, icon_url=icon_url)
        return self

    def set_channel_ctx(self, *, channel, is_channel_scope) -> Self:
        fields = []
        fields.append(f"**Channel:** {channel.mention} (`{channel.id}`)")
        fields.append(f"**Category:** {channel.category.name}")
        fields.append(f"**Was in Channel:** {is_channel_scope}")
        field = "\n".join(fields)
        self.add_field(name="📍 Channel Context", value=field, inline=True)
        return self

    def set_reason(self, reason) -> Self:
        self.add_field(name="📝 Reason", value=f"```{reason}```", inline=False)
        return self

    def set_description(self, *, channel, target) -> Self:
        if hasattr(target, "id"):
            self.description = (
                f"**Target:** {target.mention} moderated in {channel.mention}"
            )
        else:
            self.description = (
                f"**Target:** {target.get("name")} moderated in {channel.mention}"
            )
        return self


class StreamService:
    __CHUNK_SIZE = 12
    MODEL = Stream

    def __init__(
        self,
        *,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        duration_builder=None,
        emoji=None,
        moderator_service=None,
        paginator_service=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__dictionary_service = dictionary_service
        self.__database_factory.model = self.MODEL
        self.__duration_builder = duration_builder
        self.__emoji = emoji
        self.__paginator_service = paginator_service
        self.__moderator_service = moderator_service

    async def send_log(
        self,
        channel: discord.abc.GuildChannel,
        identifier: str,
        member,
        author: discord.Member | None = None,
        duration_value: str | None = None,
        is_channel_scope: bool = False,
        is_modification: bool = False,
        reason: str = "No reason provided",
        source: (
            Union[commands.Context, discord.Interaction, discord.Message] | None
        ) = None,
    ):
        message = None
        if isinstance(source, commands.Context):
            message = source.message
        elif isinstance(source, discord.Message):
            message = source
        pages = []
        executor_role = (
            await self.__moderator_service.resolve_highest_role_at_all(
                member_snowflake=int(author.id),
            )
            if author
            else "Unknown"
        )
        if hasattr(member, "id"):
            member_snowflake = member.id
        else:
            member_snowflake = member.get("id", None)
        target_role = await self.__moderator_service.resolve_highest_role_at_all(
            member_snowflake=int(member_snowflake),
        )
        embed = (
            StreamEmbed(duration_builder=self.__duration_builder)
            .set_tn(url=author.display_avatar.url)
            .set_title(identifier=identifier, is_modification=is_modification)
            .set_description(channel=channel, target=member)
            .set_target(target=member, highest_role=target_role)
            .set_executor(author=author, highest_role=executor_role)
            .set_action(duration_value=duration_value)
            .set_message_ctx(identifier=identifier, message=message)
            .set_channel_ctx(channel=channel, is_channel_scope=is_channel_scope)
            .set_reference(
                channel=channel, target=member, message=message, source=source
            )
        )
        pages.append(embed)
        embed = StreamEmbed(color=None, description=None, title=None, url=None)
        embed.set_title(
            identifier=identifier, is_modification=is_modification
        ).set_reason(reason=reason)
        pages.append(embed)

        streaming = await self.__database_factory.select(singular=False)
        for stream in streaming:
            channel_obj = self.__bot.get_channel(stream.target_channel_snowflake)
            if channel_obj:
                perms = channel_obj.permissions_for(channel_obj.guild.me)
                if perms.send_messages and not channel_obj.guild.me.is_timed_out():
                    await self.__paginator_service.start(
                        channel=channel_obj, pages=pages
                    )
        return

    async def build_dictionary(self, obj):
        streaming = []
        dictionary = {}
        if isinstance(obj, discord.Guild):
            streaming = await self.__database_factory.select(guild_snowflake=obj.id)
        elif isinstance(obj, discord.abc.GuildChannel):
            streaming = await self.__database_factory.select(channel_snowflake=obj.id)
        else:
            streaming = await self.__database_factory.select()
        if streaming:
            for stream in streaming:
                dictionary.setdefault(stream.guild_snowflake, {"channels": {}})
                dictionary[stream.guild_snowflake]["channels"].setdefault(
                    stream.target_channel_snowflake, {"sources": []}
                )
                dictionary[stream.guild_snowflake]["channels"][
                    stream.target_channel_snowflake
                ]["sources"].append(stream.source_channel_snowflake)
        return dictionary

    async def build_pages(self, is_at_home, obj):
        lines, pages = [], []

        obj_name = "All Servers"
        if obj:
            obj_name = obj.name
        title = f"{self.__emoji.get_random_emoji()} Streaming Routes for {obj_name}"

        dictionary = await self.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=StreamDictionary, dictionary=dictionary
        )

        for guild_snowflake, guild_data in processed_dictionary.data.items():
            stream_n = 0
            field_count = 0
            lines = []
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, entry in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                sources = []
                for source_snowflake in entry["sources"]:
                    source_channel = guild.get_channel(source_snowflake)
                    if source_channel:
                        sources.append(source_channel.mention)
                if sources:
                    source_lines = "\n".join(f"• {s}" for s in sources)
                else:
                    source_lines = "• All sources"
                lines.append(
                    f"**Target:** {channel.mention}\n**Sources:**\n{source_lines}"
                )
                field_count += 1
                stream_n += 1
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = self.__dictionary_service.flush_page(
                        embed, pages, title, guild.name
                    )
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            original_description = embed.description or ""
            embed.description = f"**{original_description} ({stream_n})**"
            pages.append(embed)
        if is_at_home:
            pages.extend(processed_dictionary.skipped_channels)
            pages.extend(processed_dictionary.skipped_guilds)
        return pages

    async def toggle_stream(
        self,
        source,
        target_channel,
    ):
        if source:
            stream = await self.__database_factory.select(
                singular=True,
                guild_snowflake=target_channel.guild.id,
                source_channel_snowflake=source.id,
                target_channel_snowflake=target_channel.id,
            )
            if stream:
                await self.__database_factory.delete(
                    guild_snowflake=target_channel.guild.id,
                    source_channel_snowflake=source.id,
                    target_channel_snowflake=target_channel.id,
                )
            source_text = f"from {source.mention}"
        else:
            stream = await self.__database_factory.select(
                singular=True,
                guild_snowflake=target_channel.guild.id,
                target_channel_snowflake=target_channel.id,
            )
            if stream:
                await self.__database_factory.delete(
                    guild_snowflake=target_channel.guild.id,
                    source_channel_snowflake=source.id,
                )
            source_text = "from all channels"
        if stream:
            action = "deleted"
        else:
            action = "created"
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} Tracking {action.capitalize()} {source_text} to {target_channel.mention}",
            color=0x00FF00,
        )
        return [embed]
