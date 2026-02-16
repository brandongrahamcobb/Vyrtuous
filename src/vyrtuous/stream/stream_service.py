from copy import copy

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
from typing import Self
from datetime import datetime, timezone
from typing import Union

import discord
from discord.ext import commands

from vyrtuous.stream.stream import Stream


class StreamEmbed(discord.Embed):
    def __init__(self, *, color=None, description=None, title=None, url=None):
        self.__action = None
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
            self.__action = "banned"
        elif identifier in ("flag", "unflag"):
            if is_modification:
                self.title = "🔄 User Unflagged"
            else:
                self.title = "🚩 User Flaged"
            self.__action = "flagged"
        elif identifier in ("tmute", "untmute"):
            if is_modification:
                self.title = "🔄 User Unmuted"
            else:
                self.title = "📝 User Text Muted"
            self.__action = "text muted"
        elif identifier in ("vmute", "unvmute"):
            if is_modification:
                self.title = "🔄 User Unmuted"
            else:
                self.title = "🎙️ User Voice Muted"
            self.__action = "voice muted"
        return self

    def set_thumbnail(self, *, url):
        self.set_thumbnail(url=url)

    def set_executor(self, *, author, highest_role) -> Self:
        fields = []
        fields.append(f"**Executor:** {author.display_name} (@{author.name})")
        fields.append(f"**Executor ID:** `{author.id}`")
        fields.append(f"**Top Role:** {highest_role}")
        field = "\n".join(fields)
        self.add_field(name="👮‍♂️ Executed By", value=field, inline=True)
        return self

    def set_target(self, *, target, highest_role) -> Self:
        fields = []
        fields.append(f"**Display Name:** {target.display_name}")
        fields.append(f"**Username:** @{target.name}")
        fields.append(f"**User ID:** `{target.id}")
        fields.append(f"**Account Age:** <t:{int(target.created_at.timestamp())}:R>")
        fields.append(f"**Server Join:** <t:{int(target.joined_at.timestamp())}:R>")
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

    def set_action(self, *, duration) -> Self:
        if duration is not None:
            dt = "⏱️ Temporary"
            expiration = f"{duration}"
        else:
            dt = "♾️ Permanent"
            expiration = "Never"
        fields = []
        fields.append(f"**Expires:** {expiration}")
        fields.append(f"**Alias:** {self.__action}")
        self.add_field(name=f"**Type:** {dt}", value="\n".join(fields), inline=False)
        return self

    def set_reference(self, *, channel, target, message, source) -> Self:
        text = f"Ref: {target.id}-{channel.id} | Msg: {message.id if message else 'Hidden'}"
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
        self.add_field(name=f"📝 Reason", value=f"```{reason}```", inline=False)
        return self

    def set_description(self, *, channel, target) -> Self:
        self.description = (
            f"**Target:** {target.mention} {self.action} in {channel.mention}"
        )
        return self


class StreamService:
    __CHUNK_SIZE = 7
    MODEL = Stream

    def __init__(
        self,
        *,
        author_service=None,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        duration_service=None,
        emoji=None,
        moderator_service=None,
    ):
        self.__author_service = author_service
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__dictionary_service = dictionary_service
        self.__database_factory.model = self.MODEL
        self.__emoji = emoji
        self.__duration_service = duration_service
        self.__moderator_service = moderator_service

    async def send_log(
        self,
        author: discord.Member,
        channel: discord.abc.GuildChannel,
        identifier: str,
        target: discord.Member,
        duration: str | None = None,
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
        executor_role = await self.__moderator_service.resolve_highest_role_at_all(
            member_snowflake=int(author.id),
        )
        target_role = await self.__moderator_service.resolve_highest_role_at_all(
            member_snowflake=int(target.id),
        )
        embed = (
            StreamEmbed()
            .set_title(identifier=identifier, is_modification=is_modification)
            .set_description(channel=channel, target=target)
            .set_target(target=target, highest_role=target_role)
            .set_executor(author=author, highest_role=executor_role)
            .set_action(duration=duration)
            .set_message_ctx(identifier=identifier, message=message)
            .set_channel_ctx(channel=channel, is_channel_scope=is_channel_scope)
            .set_reference(
                channel=channel, target=target, message=message, source=source
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
            channel_obj = self.__bot.get_channel(stream.channel_snowflake)
            if channel_obj:
                perms = channel_obj.permissions_for(channel_obj.guild.me)
                if perms.send_messages and not channel_obj.guild.me.is_timed_out():
                    paginator = self.__paginator_service(self.__bot, channel_obj, pages)
                    await paginator.start()
        return

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        pages = []
        dictionary = {}
        streaming = await self.__database_factory.select(singular=False, **where_kwargs)
        for stream in streaming:
            dictionary.setdefault(stream.guild_snowflake, {"channels": {}})
            dictionary[stream.guild_snowflake]["channels"][stream.channel_snowflake] = {
                "enabled": stream.enabled,
                "entry_type": stream.entry_type,
                "snowflakes": stream.snowflakes,
            }
        skipped_channels = self.__dictionary_service.generate_skipped_channels(
            dictionary
        )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        # if is_at_home:
        #     if skipped_channels:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_dict_pages(
        #                 skipped=skipped_channels,
        #                 title="Skipped Channels in Server",
        #             )
        #         )
        #     if skipped_guilds:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_set_pages(
        #                 skipped=skipped_guilds,
        #                 title="Skipped Servers",
        #             )
        #         )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Streaming Routes"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )
        embed = discord.Embed(
            title=title, description="Default view", color=discord.Color.blue()
        )
        if not dictionary:
            pages = [embed]

        stream_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, entry in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                if not channel:
                    continue
                status = "\u2705" if entry["enabled"] else "\u26d4"
                lines.append(
                    f"{status}**Channel:** {channel.mention}\n**Type:** {entry['entry_type']}"
                )
                if isinstance(
                    object_dict.get("object", None), discord.abc.GuildChannel
                ):
                    lines.append(f"**Snowflakes:** {entry['snowflakes']}")
                field_count += 1
                stream_n += 1
                if field_count >= self.__CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed = flush_page(embed, pages, title, guild.name)
                    lines = []
                    field_count = 0
            if lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)
        if pages:
            pages[0].description = f"**({stream_n})**"
        return pages

    async def modify_stream(
        self,
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
            stream = await self.__database_factory.select(singular=True, **where_kwargs)
            enabled = not stream.enabled
            action = "enabled" if enabled else "disabled"
            set_kwargs = {"enabled": enabled}
            await self.__database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )
        if action and entry_type:
            match action.lower():
                case "create" | "modify":
                    resolved_channels = []
                    if action.lower() == "create":
                        stream = self.MODEL(
                            **where_kwargs,
                            enabled=True,
                            snowflakes=resolved_channels,
                        )
                        await self.__database_factory.create(stream)
                        action = "created"
                    else:
                        set_kwargs = {"snowflakes": resolved_channels}
                        await self.__database_factory.update(
                            set_kwargs=set_kwargs, where_kwargs=where_kwargs
                        )
                        action = "modified"
                case "delete":
                    await Stream.delete(**where_kwargs)
                    action = "deleted"
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} Tracking {action.capitalize()} for {channel_dict.get('mention', None)}",
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
