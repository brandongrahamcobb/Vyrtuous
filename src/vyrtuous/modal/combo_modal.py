"""!/bin/python3
reason_modal.py The purpose of this program is to provide the reason utility modal which is used to finalize infractions.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord

from vyrtuous.aliases import (
    flag_alias_service,
    text_mute_alias_service,
    voice_mute_alias_service,
)
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.flag import Flag
from vyrtuous.db.text_mute import TextMute
from vyrtuous.db.voice_mute import VoiceMute
from vyrtuous.models.duration import DurationBuilder, DurationObject
from vyrtuous.utils.errors.error import ChannelNotFound, GuildNotFound, MemberNotFound
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.tracking import stream_service

MODELS = [Flag, TextMute, VoiceMute]

class ComboModal(discord.ui.Modal):
    def __init__(
        self,
        author_snowflake: int,
        channel_snowflake: int,
        duration: DurationObject | None,
        guild_snowflake: int,
        flag_enabled: bool,
        member_snowflake: int,
        records: list[Flag | TextMute | VoiceMute],
        text_mute_enabled: bool,
        tick: Tick,
        voice_mute_enabled: bool
    ):
        super().__init__(title="Reason", timeout=120)
        self.__author_snowflake = author_snowflake
        self.__channel_snowflake = channel_snowflake
        self.__duration = duration
        self.__flag_enabled = flag_enabled
        self.__guild_snowflake = guild_snowflake
        self.__member_snowflake = member_snowflake
        self.__records = records
        self.__text_mute_enabled = text_mute_enabled
        self.__tick = tick
        self.__voice_mute_enabled = voice_mute_enabled

    async def setup(self):
        self.reason_selection = discord.ui.TextInput(
            label="Type the reason",
            style=discord.TextStyle.paragraph,
            required=True,
            default="",
        )
        self.add_item(self.reason_selection)

    async def on_submit(self, interaction) -> None:
        if self.__duration is None:
            return
        if self.reason_selection.value == "":
            return await interaction.response.send_message(
                "A reason is required for flags."
            )
        else:
            reason = self.reason_selection.value
        self.__tick.update_source(interaction=interaction)
        await self.__tick.defer()
        duration_builder: DurationBuilder = DurationBuilder()
        flag_updated = False
        text_mute_updated = False
        voice_mute_updated = False
        for record in self.__records:
            match record:
                case Flag() as flag:
                    database_factory: DatabaseFactory = DatabaseFactory(Flag)
                    set_kwargs = {
                        "reason": reason,
                    }
                    where_kwargs = {
                        "channel_snowflake": self.__channel_snowflake,
                        "guild_snowflake": self.__guild_snowflake,
                        "member_snowflake": self.__member_snowflake,
                    }
                    await database_factory.update(
                        set_kwargs=set_kwargs, where_kwargs=where_kwargs
                    )
                    flag_updated = True
                case TextMute() as tmute:
                    expires_in = duration_builder.load(self.__duration).to_expires_in()
                    database_factory: DatabaseFactory = DatabaseFactory(TextMute)
                    set_kwargs = {
                        "expires_in": expires_in,
                        "reason": reason
                    }
                    where_kwargs = {
                        "channel_snowflake": self.__channel_snowflake,
                        "guild_snowflake": self.__guild_snowflake,
                        "member_snowflake": self.__member_snowflake,
                    }
                    await database_factory.update(
                        set_kwargs=set_kwargs, where_kwargs=where_kwargs
                    )
                    text_mute_updated = True
                case VoiceMute() as vmute:
                    target = "command"
                    expires_in = duration_builder.load(self.__duration).to_expires_in()
                    database_factory: DatabaseFactory = DatabaseFactory(VoiceMute)
                    set_kwargs = {
                        "expires_in": expires_in,
                        "reason": reason
                    }
                    where_kwargs = {
                        "channel_snowflake": self.__channel_snowflake,
                        "guild_snowflake": self.__guild_snowflake,
                        "member_snowflake": self.__member_snowflake,
                    }
                    await database_factory.update(
                        set_kwargs=set_kwargs, where_kwargs=where_kwargs
                    )
                    voice_mute_updated = True
        duration_builder: DurationBuilder = DurationBuilder()
        expires_in = (
            duration_builder.load(self.__duration).to_expires_in()
            if self.__duration
            else None
        )
        display = True
        is_channel_scope = False
        target = "command"
        for model in MODELS:
            database_factory: DatabaseFactory = DatabaseFactory(model)
            if model == Flag and not flag_updated:
                if self.__flag_enabled:
                    flag = Flag(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=reason,
                    )
                    await database_factory.create(flag)
                    is_channel_scope = await flag_alias_service.enable(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=reason,
                    )
            elif model == TextMute and not text_mute_updated:
                if self.__text_mute_enabled:
                    tmute = TextMute(
                        channel_snowflake=self.__channel_snowflake,
                        expires_in=expires_in,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=reason,
                    )
                    await database_factory.create(tmute)
                    is_channel_scope = await text_mute_alias_service.enable(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=reason,
                    )
            elif model == VoiceMute and not voice_mute_updated:
                if self.__voice_mute_enabled:
                    vmute = VoiceMute(
                        channel_snowflake=self.__channel_snowflake,
                        expires_in=expires_in,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=reason,
                        target=target,
                    )
                    await database_factory.create(vmute)
                    is_channel_scope = await voice_mute_alias_service.enable(
                        channel_snowflake=self.__channel_snowflake,
                        guild_snowflake=self.__guild_snowflake,
                        member_snowflake=self.__member_snowflake,
                        reason=reason,
                    )
        await stream_service.log(
            author_snowflake=self.__author_snowflake,
            channel_snowflake=self.__channel_snowflake,
            display=display,
            duration=self.__duration,
            guild_snowflake=self.__guild_snowflake,
            identifier="combo",
            is_channel_scope=is_channel_scope,
            member_snowflake=self.__member_snowflake,
            message_snowflake=None,
            message_channel_snowflake=None,
            reason=reason,
            role_snowflake=None,
            target=target,
        )
        embed = build_combo_embed(
            channel_snowflake=self.__channel_snowflake,
            duration=self.__duration,
            flag_enabled=self.__flag_enabled,
            flag_updated=flag_updated,
            guild_snowflake=self.__guild_snowflake,
            member_snowflake=self.__member_snowflake,
            reason=reason,
            text_mute_enabled=self.__text_mute_enabled,
            text_mute_updated=text_mute_updated,
            voice_mute_enabled=self.__voice_mute_enabled,
            voice_mute_updated=voice_mute_updated,
        )
        await self.__tick.end(success=embed)

    async def on_error(
        self,
        interaction: discord.Interaction,
        error: Exception,
    ):
        await self.__tick.end(error=str(error), ephemeral=True)


def build_combo_embed(
    channel_snowflake: int,
    duration: DurationObject,
    flag_enabled: bool,
    flag_updated: bool,
    guild_snowflake: int,
    member_snowflake: int,
    reason: str,
    text_mute_enabled: bool,
    text_mute_updated: bool,
    voice_mute_enabled: bool,
    voice_mute_updated: bool,
) -> discord.Embed:
    fields = []
    if flag_enabled:
        fields.append("flagged")
    if text_mute_enabled:
        fields.append("text-muted")
    if voice_mute_enabled:
        fields.append("voice-muted")
    bot: DiscordBot = DiscordBot.get_instance()
    duration_builder = DurationBuilder()
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise ChannelNotFound(str(channel_snowflake))
    member = guild.get_member(member_snowflake)
    if member:
        display_name = member.display_name
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(member_snowflake)
        if not simplified_member:
            raise MemberNotFound(str(member_snowflake))
        display_name = simplified_member[0]
        member_str = display_name
    embed = discord.Embed(
        title=f"{emojis.get_random_emoji()} {display_name} was {", ".join(fields)}",
        description=(
            f"**User:** {member_str}\n"
            f"**Channel:** {channel.mention}\n"
            f"**Expires:** {duration_builder.load(duration=duration).to_unix_ts()}\n"
            f"**Reason:** {reason}"
        ),
        color=discord.Color.blue(),
    )
    if member:
        embed.set_thumbnail(url=member.display_avatar.url)
    return embed
