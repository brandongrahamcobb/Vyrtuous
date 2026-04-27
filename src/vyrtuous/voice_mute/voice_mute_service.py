"""!/bin/python3
voice_mute_service.py The purpose of this program is to extend AliasService to service the voice mute infraction.

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
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Union

import discord
from discord.ext import commands

from vyrtuous.active_members import active_member_service
from vyrtuous.moderator.moderator_service import HasEqualOrLowerRole
from vyrtuous.voice_mute.voice_mute import VoiceMute


@dataclass
class VoiceMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class VoiceMuteService:
    __CHUNK_SIZE = 12
    MODEL = VoiceMute
    voice_muted_members = {}

    def __init__(
        self,
        *,
        active_member_service=None,
        bot=None,
        database_factory=None,
        data_service=None,
        dictionary_service=None,
        duration_builder=None,
        emoji=None,
        moderator_service=None,
        stream_service=None,
    ):
        self.__active_member_service = active_member_service
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__data_service = data_service
        self.__dictionary_service = dictionary_service
        self.__duration_builder = duration_builder
        self.__emoji = emoji
        self.__moderator_service = moderator_service
        self.__stream_service = stream_service

    async def populate(self):
        voice_muted_members = await self.__database_factory.select()
        for voice_muted_member in voice_muted_members:
            guild = self.__bot.get_guild(voice_muted_member.guild_snowflake)
            if not guild:
                continue
            self.text_muted_members[voice_muted_member.member_snowflake] = {
                "last_active": None,
                "name": voice_muted_member.display_name,
            }

    async def enforce_or_undo(
        self,
        ctx,
        default_ctx,
        source: Union[commands.Context, discord.Interaction, discord.Message],
        state,
    ):
        obj = await self.__database_factory.select(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member.id,
            singular=True,
        )
        if obj:
            await self.undo(
                ctx=ctx, default_ctx=default_ctx, source=source, state=state
            )
        else:
            await self.enforce(
                ctx=ctx, default_ctx=default_ctx, source=source, state=state
            )

    async def clean_expired(self):
        expired_voice_mutes = await self.__database_factory.select(expired=True)
        if expired_voice_mutes:
            for expired_voice_mute in expired_voice_mutes:
                channel_snowflake = int(expired_voice_mute.channel_snowflake)
                guild_snowflake = int(expired_voice_mute.guild_snowflake)
                member_snowflake = int(expired_voice_mute.member_snowflake)
                target = expired_voice_mute.target
                guild = self.__bot.get_guild(guild_snowflake)
                kwargs = {
                    "channel_snowflake": channel_snowflake,
                    "guild_snowflake": guild_snowflake,
                    "member_snowflake": member_snowflake,
                    "target": target,
                }
                if guild is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired voice-mute."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired voice-mute."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.name}), cleaning up expired voice-mute."
                    )
                    continue
                await self.__database_factory.delete(**kwargs)
                if (
                    member.voice
                    and member.voice.channel
                    and member.voice.channel.id == channel_snowflake
                ):
                    try:
                        await member.edit(mute=False)
                        self.__bot.logger.info(
                            f"Undone voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake})."
                        )
                    except discord.Forbidden as e:
                        self.__bot.logger.warning(
                            f"Unable to undo voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}). {str(e).capitalize()}"
                        )
                else:
                    self.__bot.logger.info(
                        f"Member {member.display_name} ({member.id}) is not in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), skipping undo voice-mute."
                    )

    async def build_dictionary(self, obj):
        voice_mutes = []
        dictionary = {}
        if isinstance(obj, discord.Guild):
            voice_mutes = await self.__database_factory.select(guild_snowflake=obj.id)
        elif isinstance(obj, discord.abc.GuildChannel):
            voice_mutes = await self.__database_factory.select(channel_snowflake=obj.id)
        elif isinstance(obj, discord.Member):
            voice_mutes = await self.__database_factory.select(member_snowflake=obj.id)
        else:
            voice_mutes = await self.__database_factory.select()
        if voice_mutes:
            for voice_mute in voice_mutes:
                dictionary.setdefault(voice_mute.guild_snowflake, {"members": {}})
                dictionary[voice_mute.guild_snowflake]["members"].setdefault(
                    voice_mute.member_snowflake, {"voice_mutes": {}}
                )
                dictionary[voice_mute.guild_snowflake]["members"][
                    voice_mute.member_snowflake
                ]["voice_mutes"].setdefault(voice_mute.channel_snowflake, {})
                dictionary[voice_mute.guild_snowflake]["members"][
                    voice_mute.member_snowflake
                ]["voice_mutes"][voice_mute.channel_snowflake].update(
                    {"reason": voice_mute.reason, "expires_in": voice_mute.expires_in}
                )
        return dictionary

    async def build_pages(self, is_at_home, obj):
        lines, pages = [], []

        obj_name = "All Servers"
        if not isinstance(obj, int):
            obj_name = obj.name
        else:
            member = self.__active_member_service.active_member.get(obj, None)
            if member:
                obj_name = member.get("name", None)
            else:
                return "No active voice-mutes found."
        title = f"{self.__emoji.get_random_emoji()} Voice Mutes for {obj_name}"

        dictionary = await self.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=VoiceMuteDictionary, dictionary=dictionary
        )

        vmute_n = 0
        for guild_snowflake, guild_data in processed_dictionary.data.items():
            field_count = 0
            lines = []
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, voice_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if member:
                    if not isinstance(obj, discord.Member):
                        lines.append(
                            f"**User:** {member.display_name} {member.mention}"
                        )
                        field_count += 1
                    elif not thumbnail:
                        embed.set_thumbnail(url=obj.display_avatar.url)
                        thumbnail = True
                else:
                    member = self.__active_member_service.active_members.get(
                        member_snowflake, None
                    )
                    if member:
                        display_name = member.get("name", None)
                        lines.append(f"**User:** {display_name} ({member_snowflake})")
                for channel_snowflake, channel_dictionary in voice_mute_dictionary.get(
                    "voice_mutes", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(obj, discord.abc.GuildChannel):
                        lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(obj, discord.Member):
                        lines.append(
                            f"**Expires in:** {channel_dictionary['expires_in']}"
                        )
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
                    vmute_n += 1
                    field_count += 1
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
            pages.append(embed)
        if pages:
            original_description = embed.description or ""
            embed.description = f"**{original_description} ({vmute_n})**"
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)
        if not pages:
            return "No voice mutes found."
        return pages

    async def room_mute(self, author, channel, reason):
        muted_members, pages, skipped_members, failed_members = [], [], [], []
        for member in channel.members:
            if member.id == author.id:
                continue
            voice_mute = await self.__database_factory.select(
                channel_snowflake=channel.id,
                member_snowflake=member.id,
                target="user",
                singular=True,
            )
            try:
                await self.__moderator_service.has_equal_or_lower_role(
                    channel_snowflake=channel.id,
                    guild_snowflake=channel.guild.id,
                    member_snowflake=author.id,
                    target_member_snowflake=member.id,
                )
            except:
                skipped_members.append(member)
                continue
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel.id:
                    try:
                        await member.edit(mute=True)
                    except Exception as e:
                        self.__bot.logger.warning(
                            f"Unable to voice-mute member "
                            f"{member.display_name} ({member.id}) in channel "
                            f"{channel.name} ({channel.id}) in guild "
                            f"{channel.guild.name} ({channel.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
            voice_mute = self.MODEL(
                channel_snowflake=channel.id,
                display_name=member.display_name,
                expires_in=expires_in,
                guild_snowflake=channel.guild.id,
                member_snowflake=member.id,
                reason=reason,
                target="user",
            )
            await self.__database_factory.create(voice_mute)
            muted_members.append(member)
        description_lines = [
            f"**Channel:** {channel.mention}",
            f"**Muted:** {len(muted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f"**Skipped:** {
                len(channel.members) - len(muted_members) - len(failed_members)
            }",
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{self.__emoji.get_random_emoji()} Room Mute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)
        return pages

    async def room_unmute(self, channel):
        unmuted_members, pages, skipped_members, failed_members = [], [], [], []
        for member in channel.members:
            voice_mute = await self.__database_factory.select(
                target="user",
                channel_snowflake=channel.id,
                member_snowflake=member.id,
                singular=True,
            )
            if not voice_mute:
                skipped_members.append(member)
                continue
            await self.__database_factory.delete(
                target="user", channel_snowflake=channel.id, member_snowflake=member.id
            )
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel.id:
                    try:
                        await member.edit(mute=False)
                    except Exception as e:
                        self.__bot.logger.warning(
                            f"Unable to undo voice-mute "
                            f"for member {member.display_name} ({member.id}) "
                            f"in channel {channel.name} ({channel.id}) "
                            f"in guild {channel.guild.name} "
                            f"({channel.guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            unmuted_members.append(member)
        description_lines = [
            f"**Channel:** {channel.mention}",
            f"**Unmuted:** {len(unmuted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f"**Skipped:** {
                len(channel.members) - len(unmuted_members) - len(failed_members)
            }",
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{self.__emoji.get_random_emoji()} Room Unmute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)
        return pages

    async def delete(
        self,
        guild_snowflake=None,
        channel_snowflake=None,
        member_snowflake=None,
    ):
        kwargs = {}
        if channel_snowflake:
            kwargs.update({"channel_snowflake": channel_snowflake})
        if guild_snowflake:
            kwargs.update({"guild_snowflake": guild_snowflake})
        if member_snowflake:
            kwargs.update({"member_snowflake": member_snowflake})
        objects = await self.__database_factory.select(**kwargs)
        for obj in objects:
            await self.__database_factory.delete_by_cls(obj, **kwargs)
            guild = self.__bot.get_guild(obj.guild_snowflake)
            channel = guild.get_channel(obj.channel_snowflake)
            member = guild.get_member(obj.member_snowflake)
            is_channel_scope = False
            if member.voice and member.voice.channel:
                is_channel_scope = True
            await self.__stream_service.send_log(
                channel=channel,
                identifier="unvmute",
                is_channel_scope=is_channel_scope,
                member=member,
                reason="Right-click unvoice-mute.",
            )
            await self.__data_service.save_data(
                channel=channel,
                identifier="unvmute",
                member=member,
            )

    async def create(
        self,
        guild_snowflake=None,
        channel_snowflake=None,
        member_snowflake=None,
    ):
        kwargs = {}
        if channel_snowflake:
            kwargs.update({"channel_snowflake": channel_snowflake})
        if guild_snowflake:
            kwargs.update({"guild_snowflake": guild_snowflake})
        if member_snowflake:
            kwargs.update({"member_snowflake": member_snowflake})
        objects = await self.__database_factory.select(**kwargs)
        for obj in objects:
            await self.__database_factory.creat(obj, **kwargs)
            guild = self.__bot.get_guild(obj.guild_snowflake)
            channel = guild.get_channel(obj.channel_snowflake)
            member = guild.get_member(obj.member_snowflake)
            is_channel_scope = False
            if member.voice and member.voice.channel:
                is_channel_scope = True
            await self.__stream_service.send_log(
                channel=channel,
                identifier="vmute",
                is_channel_scope=is_channel_scope,
                member=member,
                reason="Right-click voice-mute.",
            )
            await self.__data_service.save_data(
                channel=channel,
                identifier="vmute",
                member=member,
            )

    async def enforce_log(
        self, author, channel, duration_value, is_channel_scope, member, source, reason
    ):
        await self.__stream_service.send_log(
            author=author,
            channel=channel,
            duration_value=duration_value,
            identifier="vmute",
            is_channel_scope=is_channel_scope,
            member=member,
            source=source,
            reason=reason,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="vmute",
            duration_value=duration_value,
            reason=reason,
            member=member,
        )

    async def enforce(self, ctx, default_ctx, source, state):
        voice_mute = self.MODEL(
            channel_snowflake=ctx.channel.id,
            display_name=ctx.display_name,
            expires_in=ctx.expires_in,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member_snowflake,
            reason=ctx.reason,
            target="user",
        )
        await self.__database_factory.create(voice_mute)
        is_channel_scope = False
        member = ctx.guild.get_member(ctx.member_snowflake)
        if ctx.member.voice and ctx.member.voice.channel and member:
            if ctx.member.voice.channel.id == ctx.channel.id:
                is_channel_scope = True
                try:
                    await ctx.member.edit(mute=True, reason=ctx.reason)
                except discord.Forbidden as e:
                    self.__bot.logger.info(str(e).capitalize())
                    return await state.end(error=str(e).capitalize())
        if not member:
            member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
        self.voice_muted_members.update(
            {ctx.member_snowflake: {"name": ctx.display_name}}
        )
        await self.enforce_log(
            author=default_ctx.author,
            channel=ctx.channel,
            duration_value=ctx.duration_value,
            is_channel_scope=is_channel_scope,
            member=member,
            source=source,
            reason=ctx.reason,
        )
        embed = await self.act_embed(ctx=ctx)
        return await state.end(success=embed)

    async def undo_log(self, author, channel, is_channel_scope, member, source):
        await self.__stream_service.send_log(
            author=author,
            channel=channel,
            identifier="unvmute",
            is_channel_scope=is_channel_scope,
            is_modification=True,
            member=member,
            source=source,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="unvmute",
            is_modification=True,
            member=member,
        )

    async def undo(self, ctx, default_ctx, source, state):
        await self.__database_factory.delete(
            channel_snowflake=ctx.channel.id,
            display_name=ctx.display_name,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member_snowflake,
        )
        member = ctx.guild.get_member(ctx.member_snowflake)
        is_channel_scope = False
        if ctx.member.voice and ctx.member.voice.channel and member:
            try:
                is_channel_scope = True
                await ctx.member.edit(mute=False)
            except discord.Forbidden as e:
                self.__bot.logger.info(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        if not member:
            member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
        # del self.voice_muted_members[ctx.member_snowflake]
        await self.undo_log(
            author=default_ctx.author,
            channel=ctx.channel,
            is_channel_scope=is_channel_scope,
            member=member,
            source=source,
        )
        embed = await self.undo_embed(ctx=ctx)
        return await state.end(success=embed)

    async def act_embed(self, ctx):
        member = ctx.guild.get_member(ctx.member_snowflake)
        if member:
            member_display_name = member.display_name
            member_str = member.mention
        else:
            simplified_member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
            member_display_name = simplified_member.get("name", None)
            member_str = simplified_member.get("name", None)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member_display_name} has been voice-muted",
            description=(
                f"**User:** {member_str}\n"
                f"**Channel:** {ctx.channel.mention}\n"
                f"**Expires:** {self.__duration_builder.parse(ctx.duration_value).to_unix_ts()}\n"
                f"**Reason:** {ctx.reason}"
            ),
            color=discord.Color.blue(),
        )
        if member:
            embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def undo_embed(self, ctx):
        member = ctx.guild.get_member(ctx.member_snowflake)
        if member:
            member_display_name = member.display_name
            member_str = member.mention
        else:
            simplified_member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
            member_display_name = simplified_member.get("name", None)
            member_str = simplified_member.get("name", None)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member_display_name} has been unmuted",
            description=(f"**User:** {member_str}\n**Channel:** {ctx.channel.mention}"),
            color=discord.Color.yellow(),
        )
        if member:
            embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def clean_expired_stage(self, channel, guild):
        voice_mutes = await self.__database_factory.select(
            channel_snowflake=channel.id,
            guild_snowflake=guild.id,
            target="room",
        )
        for voice_mute in voice_mutes:
            member_snowflake = voice_mute.member_snowflake
            member = guild.get_member(member_snowflake)
            if member is None:
                await self.__database_factory.delete(
                    channel_snowflake=channel.id,
                    member_snowflake=int(member_snowflake),
                    guild_snowflake=guild.id,
                    target="room",
                )
                self.__bot.logger.info(
                    f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.id}) from expired stage."
                )
                continue
            await self.__database_factory.delete(
                channel_snowflake=channel.id,
                member_snowflake=member.id,
                guild_snowflake=guild.id,
                target="room",
            )
            if (
                member.voice
                and member.voice.channel
                and member.voice.mute
                and member.voice_channel.id == channel.id
            ):
                try:
                    await member.edit(
                        mute=False, reason="Stage room closed automatically."
                    )
                    self.__bot.logger.info(
                        f"Undone voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in in guild {guild.name} ({guild.id}) after stage expired."
                    )
                except discord.Forbidden as e:
                    self.__bot.logger.warning(
                        f"Unable to undo voice-mute for member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.id}). {str(e).capitalize()}"
                    )
            else:
                self.__bot.logger.info(
                    f"Member {member.display_name} ({member.id}) is not in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild.id}), skipping undo voice-mute."
                )

    async def off_stage(self, channel):
        failed, succeeded = [], []
        for member in channel.members:
            await self.__database_factory.delete(
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                member_snowflake=member.id,
                target="room",
            )
            voice_mute = await self.__database_factory.select(
                channel_snowflake=channel.id,
                guild_snowflake=channel.guild.id,
                member_snowflake=member.id,
                target="user",
                singular=True,
            )
            if not voice_mute and member.voice and member.voice.mute:
                try:
                    await member.edit(
                        mute=False,
                        reason="Stage ended — no user-specific mute found",
                    )
                    succeeded.append(member)
                except discord.Forbidden as e:
                    self.__bot.logger.warning(
                        f"Unable to undo voice-mute "
                        f"for member {member.display_name} ({member.id}) in "
                        f"channel {channel.name} ({channel.id}) in "
                        f"guild {channel.guild.name} ({channel.guild.id}). "
                        f"{str(e).capitalize()}"
                    )
                    failed.append(member)
        return failed, succeeded

    async def on_stage(self, channel, context, duration_value):
        failed, skipped, succeeded = [], [], []
        for member in channel.members:
            if (
                await self.__moderator_service.check(
                    channel_snowflake=channel.id,
                    guild_snowflake=channel.guild.id,
                    lowest_role="Coordinator",
                )
                or member.id == context.author.id
            ):
                skipped.append(member)
                continue
            voice_mute = self.MODEL(
                channel_snowflake=channel.id,
                display_name=member.display_name,
                expires_in=self.__duration_builder.parse(
                    value=duration_value
                ).to_expires_in(),
                guild_snowflake=channel.guild.id,
                member_snowflake=member.id,
                target="room",
                reason="Stage mute",
            )
            await self.__database_factory.create(voice_mute)
            try:
                if member.voice and member.voice.channel.id == channel.id:
                    await member.edit(mute=True)
                succeeded.append(member)
            except Exception as e:
                self.__bot.logger.warning(
                    f"Unable to voice-mute "
                    f"member {member.display_name} ({member.id}) "
                    f"in channel {channel.name} ({channel.id}) "
                    f"in guild {channel.guild.name} ({channel.guild.id}). "
                    f"{str(e).capitalize()}"
                )
                failed.append(member)
        return failed, skipped, succeeded

    async def migrate(self, kwargs):
        self.__database_factory.update(**kwargs)

    async def is_voice_muted(self, channel, member):
        voice_mute = await self.__database_factory.select(
            channel_snowflake=channel.id, member_snowflake=member.id
        )
        if voice_mute:
            return True
        return False

    async def mute(self, channel, duration_value, member, target):
        expires_in = self.__duration_builder.parse(duration_value).to_expires_in()
        voice_mute = self.MODEL(
            channel_snowflake=channel.id,
            display_name=member.display_name,
            expires_in=expires_in,
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
            reason="No reason provided.",
            target=target,
        )
        await self.__database_factory.create(obj=voice_mute)
        await self.__stream_service.send_log(
            channel=channel,
            duration_value=duration_value,
            identifier="vmute",
            is_channel_scope=True,
            member=member,
            reason="Right-click voice-mute.",
        )
        await self.__data_service.save_data(
            channel=channel,
            identifier="vmute",
            duration_value=duration_value,
            member=member,
        )
        await self.toggle_mute(channel=channel, member=member, should_be_muted=True)

    async def unmute(self, channel, member, target):
        await self.__database_factory.delete(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
            target=target,
        )
        await self.__stream_service.send_log(
            channel=channel,
            identifier="unvmute",
            is_channel_scope=True,
            member=member,
            reason="Right-click unvoice-mute.",
        )
        await self.__data_service.save_data(
            channel=channel,
            identifier="unvmute",
            member=member,
        )
        await self.toggle_mute(channel=channel, member=member, should_be_muted=False)

    async def toggle_mute(self, channel, member, should_be_muted):
        try:
            await member.edit(
                mute=should_be_muted,
                reason=f"Setting mute to {should_be_muted} in {channel.name}",
            )
        except discord.Forbidden as e:
            self.__bot.logger.warning(
                f"No permission to edit "
                f"mute for {member.display_name}. {str(e).capitalize()}"
            )
        except discord.HTTPException as e:
            self.__bot.logger.warning(
                f"Failed to edit mute for {member.display_name}: {str(e).capitalize()}"
            )

    async def match(self):
        objects = await self.__database_factory.select()
        object_lookup = {
            (obj.guild_snowflake, obj.channel_snowflake, obj.member_snowflake)
            for obj in objects
        }
        async with self.__bot.db_pool.acquire() as conn:
            logs = await conn.fetch("""
                SELECT DISTINCT ON (guild_snowflake, channel_snowflake, target_snowflake)
                    guild_snowflake,
                    channel_snowflake,
                    target_snowflake,
                    identifier,
                    created_at
                FROM moderation_logs
                WHERE identifier IN ('vmute', 'unvmute')
                ORDER BY
                    guild_snowflake,
                    channel_snowflake,
                    target_snowflake,
                    created_at DESC
            """)
        for log in logs:
            key = (
                log["guild_snowflake"],
                log["channel_snowflake"],
                log["target_snowflake"],
            )
            if log["identifier"] == "vmute" and key not in object_lookup:
                guild = self.__bot.get_guild(log["guild_snowflake"])
                if not guild:
                    continue
                channel = guild.get_channel(log["channel_snowflake"])
                member = guild.get_member(log["target_snowflake"])
                if channel and member:
                    await self.__data_service.save_data(
                        channel=channel,
                        identifier="unvmute",
                        member=member,
                    )
