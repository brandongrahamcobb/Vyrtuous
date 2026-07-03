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

from datetime import datetime, timezone
from typing import Self

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.models.duration import DurationBuilder


class StreamEmbed(discord.Embed):
    def __init__(
        self,
        *,
        color=None,
        description=None,
        title=None,
        url=None,
    ):
        self.__duration_builder = DurationBuilder()
        self.__timestamp = datetime.now(timezone.utc)
        super().__init__(
            color=color,
            description=description,
            timestamp=self.__timestamp,
            title=title,
            url=url,
        )

    def set_role(self, *, guild_snowflake: int, role_snowflake: int) -> Self:
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        role = guild.get_role(role_snowflake)
        if role is None:
            raise commands.RoleNotFound(str(role_snowflake))
        self.add_field(name="Role", value=role.mention, inline=False)
        return self

    def set_title(self, *, identifier: str):
        self.color = 0xDC143C
        match identifier:
            case "ban":
                self.title = "🔨 User Banned"
            case "unban":
                self.title = "🔄 User Unbanned"
            case "flag":
                self.title = "🚩 User Flaged"
            case "unflag":
                self.title = "🔄 User Unflagged"
            case "tmute":
                self.title = "📝 User Text Muted"
            case "untmute":
                self.title = "🔄 User Unmuted"
            case "vmute":
                self.title = "🎙️ User Voice Muted"
            case "unvmute":
                self.title = "🔄 User Unmuted"
        return self

    def set_tn(self, *, url):
        self.set_thumbnail(url=url)
        return self

    def set_executor(
        self, *, author_snowflake: int, guild_snowflake: int, highest_role: str
    ) -> Self:
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        author = guild.get_member(author_snowflake)
        if author is None:
            raise commands.MemberNotFound(str(author_snowflake))
        fields = []
        if author:
            fields.append(f"**Executor:** {author.display_name} (@{author.name})")
            fields.append(f"**Executor id:** `{author.id}`")
        else:
            fields.append("**Executor:** Unknown")
            fields.append("**Executor id:** Unknown")
        fields.append(f"**Top Role:** {highest_role}")
        field = "\n".join(fields)
        self.add_field(name="👮‍♂️ Executed By", value=field, inline=True)
        return self

    def set_target(
        self,
        *,
        guild_snowflake: int,
        target: str | None,
        target_snowflake: int,
        highest_role: str,
    ) -> Self:
        fields = []
        if target:
            fields.append(f"Type: {target}")
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        member = guild.get_member(target_snowflake)
        if member is None:
            simplified_member = bot.registry.get(MemberState).active.get(
                target_snowflake, None
            )
            if simplified_member:
                fields.append(f"**Display Name:** {simplified_member[0]}")
                fields.append(f"**User ID:** `{target_snowflake}`")
            else:
                raise commands.MemberNotFound(str(target_snowflake))
        else:
            fields.append(f"**Display Name:** {member.display_name}")
            fields.append(f"**Username:** @{member.name}")
            fields.append(f"**User ID:** `{member.id}`")
            fields.append(
                f"**Account Age:** <t:{int(member.created_at.timestamp())}:R>"
            )
            if member.joined_at:
                fields.append(
                    f"**Server Join:** <t:{int(member.joined_at.timestamp())}:R>"
                )
        fields.append(f"**Top Role:** {highest_role}")
        field = "\n".join(fields)
        self.add_field(name="👤 Target User", value=field, inline=False)
        return self

    async def set_message_ctx(
        self,
        *,
        guild_snowflake: int,
        identifier: str,
        message_snowflake: int,
        message_channel_snowflake: int,
    ) -> Self:
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if not guild:
            raise commands.GuildNotFound(str(guild_snowflake))
        channel = guild.get_channel(message_channel_snowflake)
        if not channel or not isinstance(
            channel, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
        ):
            raise commands.ChannelNotFound(str(message_channel_snowflake))
        try:
            message = await channel.fetch_message(message_snowflake)
        except (discord.NotFound, discord.Forbidden, discord.HTTPException):
            raise
        fields = []
        fields.append(f"**Original Message ID:** `{message.id}`")
        fields.append(f"**Message Link:** [Jump to Message]({message.jump_url})")
        fields.append(f"**Command Channel:** {channel.mention}")
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

    def set_reference(
        self,
        *,
        channel_snowflake: int,
        guild_snowflake: int,
        target_snowflake: int,
        message_snowflake: int | None,
    ) -> Self:
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        text = f"Ref: {target_snowflake}-{channel_snowflake} | Msg: {message_snowflake if message_snowflake else 'Hidden'}"
        if guild.icon:
            icon_url = guild.icon.url
            self.set_footer(text=text, icon_url=icon_url)
        else:
            self.set_footer(text=text)
        return self

    def set_channel_ctx(
        self, *, channel_snowflake: int, guild_snowflake: int, is_channel_scope: bool
    ) -> Self:
        fields = []
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        channel = bot.get_channel(channel_snowflake)
        if channel is None or not isinstance(
            channel, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
        ):
            raise commands.ChannelNotFound(str(channel_snowflake))
        fields.append(f"**Channel:** {channel.mention} (`{channel.id}`)")
        if channel.category:
            fields.append(f"**Category:** {channel.category.name}")
        fields.append(f"**Was in Channel:** {is_channel_scope}")
        field = "\n".join(fields)
        self.add_field(name="📍 Channel Context", value=field, inline=True)
        return self

    def set_reason(self, reason) -> Self:
        self.add_field(name="📝 Reason", value=f"```{reason}```", inline=False)
        return self

    def set_description(
        self, *, channel_snowflake: int, guild_snowflake: int, target_snowflake: int
    ) -> Self:
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise commands.GuildNotFound(str(guild_snowflake))
        channel = bot.get_channel(channel_snowflake)
        if channel is None or not isinstance(
            channel, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
        ):
            raise commands.ChannelNotFound(str(channel_snowflake))
        member = guild.get_member(target_snowflake)
        if member is None:
            simplified_member = bot.registry.get(MemberState).active.get(
                target_snowflake, None
            )
            if simplified_member:
                display_name = simplified_member[0]
                self.description = (
                    f"**Target:** {display_name} moderated in {channel.mention}"
                )
            else:
                raise commands.MemberNotFound(str(target_snowflake))

        else:
            self.description = (
                f"**Target:** {member.mention} moderated in {channel.mention}"
            )
        return self
