"""!/bin/python3
text_mute_service.py The purpose of this program is to extend AliasService to service the text mute infraction.

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

from vyrtuous.text_mute.text_mute import TextMute


@dataclass
class TextMuteDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class TextMuteService:
    __CHUNK_SIZE = 12
    MODEL = TextMute
    text_muted_members = {}

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
        self.__stream_service = stream_service

    async def populate(self):
        text_muted_members = await self.__database_factory.select()
        for text_muted_member in text_muted_members:
            guild = self.__bot.get_guild(text_muted_member.guild_snowflake)
            if not guild:
                continue
            self.text_muted_members[text_muted_member.member_snowflake] = {
                "last_active": None,
                "name": text_muted_member.display_name,
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
        expired_text_mutes = await self.__database_factory.select(expired=True)
        if expired_text_mutes:
            for expired_text_mute in expired_text_mutes:
                channel_snowflake = int(expired_text_mute.channel_snowflake)
                guild_snowflake = int(expired_text_mute.guild_snowflake)
                member_snowflake = int(expired_text_mute.member_snowflake)
                kwargs = {
                    "channel_snowflake": channel_snowflake,
                    "guild_snowflake": guild_snowflake,
                    "member_snowflake": member_snowflake,
                }
                guild = self.__bot.get_guild(guild_snowflake)
                if guild is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired text-mute."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}), cleaning up expired text-mute."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), cleaning up expired text-mute."
                    )
                    continue
                await self.__database_factory.delete(**kwargs)
                try:
                    await channel.set_permissions(
                        target=member,
                        send_messages=None,
                        add_reactions=None,
                        reason="Cleaning up expired text-mute",
                    )
                except discord.Forbidden as e:
                    self.__bot.logger.error(str(e).capitalize())

    async def clean_overwrites(self):
        now = datetime.now(timezone.utc)
        text_mutes = await self.__database_factory.select()
        for text_mute in text_mutes:
            channel_snowflake = int(text_mute.channel_snowflake)
            guild_snowflake = int(text_mute.guild_snowflake)
            member_snowflake = int(text_mute.member_snowflake)
            where_kwargs = {
                "channel_snowflake": channel_snowflake,
                "guild_snowflake": guild_snowflake,
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"reset": True}
            if not text_mute.reset and text_mute.last_muted < now - timedelta(weeks=1):
                guild = self.__bot.get_guild(guild_snowflake)
                if guild is None:
                    self.__bot.logger.info(
                        f"Unable to locate guild {guild_snowflake} for removing overwrite."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    self.__bot.logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    self.__bot.logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}) for removing overwrite."
                    )
                    continue
                try:
                    await channel.set_permissions(
                        target=member,
                        overwrite=None,
                        reason="Resetting text-mute overwrite",
                    )
                except discord.Forbidden as e:
                    self.__bot.logger.error(str(e).capitalize())
                await self.__database_factory.update(
                    set_kwargs=set_kwargs, where_kwargs=where_kwargs
                )

    async def build_dictionary(self, obj):
        text_mutes = []
        dictionary = {}
        if isinstance(obj, discord.Guild):
            text_mutes = await self.__database_factory.select(guild_snowflake=obj.id)
        elif isinstance(obj, discord.abc.GuildChannel):
            text_mutes = await self.__database_factory.select(channel_snowflake=obj.id)
        elif isinstance(obj, discord.Member):
            text_mutes = await self.__database_factory.select(member_snowflake=obj.id)
        else:
            text_mutes = await self.__database_factory.select()
        if text_mutes:
            for text_mute in text_mutes:
                dictionary.setdefault(text_mute.guild_snowflake, {"members": {}})
                dictionary[text_mute.guild_snowflake]["members"].setdefault(
                    text_mute.member_snowflake, {"text_mutes": {}}
                )
                dictionary[text_mute.guild_snowflake]["members"][
                    text_mute.member_snowflake
                ]["text_mutes"].setdefault(text_mute.channel_snowflake, {})
                dictionary[text_mute.guild_snowflake]["members"][
                    text_mute.member_snowflake
                ]["text_mutes"][text_mute.channel_snowflake].update(
                    {"reason": text_mute.reason, "expires_in": text_mute.expires_in}
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
                return "No active text-mutes found."

        title = f"{self.__emoji.get_random_emoji()} Text Mutes for {obj_name}"

        dictionary = await self.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=TextMuteDictionary, dictionary=dictionary
        )

        for guild_snowflake, guild_data in processed_dictionary.data.items():
            tmute_n = 0
            field_count = 0
            lines = []
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, text_mute_dictionary in guild_data.get(
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

                for channel_snowflake, channel_dictionary in text_mute_dictionary.get(
                    "text_mutes", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(obj, discord.abc.GuildChannel):
                        lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(obj, discord.Member):
                        lines.append(
                            f"**Expires in:** {channel_dictionary['expires_in']}"
                        )
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
                    tmute_n += 1
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
            original_description = embed.description or ""
            embed.description = f"**{original_description} ({tmute_n})**"
            pages.append(embed)
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)
        return pages

    async def delete(
        self,
        author,
        source,
        *,
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
            if not member:
                continue
            if channel:
                try:
                    await channel.set_permissions(
                        target=member,
                        send_messages=None,
                        add_reactions=None,
                        reason="Clear command",
                    )
                except discord.Forbidden as e:
                    self.__bot.logger.error(str(e).capitalize())
            else:
                continue
            await self.undo_log(
                author=author,
                channel=channel,
                member=member,
                source=source,
            )

    async def enforce_log(
        self, author, channel, duration_value, member, source, reason
    ):
        await self.__stream_service.send_log(
            author=author,
            channel=channel,
            duration_value=duration_value,
            identifier="tmute",
            member=member,
            source=source,
            reason=reason,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="tmute",
            duration_value=duration_value,
            reason=reason,
            member=member,
        )

    async def enforce(self, ctx, default_ctx, source, state):
        text_mute = self.MODEL(
            channel_snowflake=ctx.channel.id,
            display_name=ctx.display_name,
            expires_in=ctx.expires_in,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member_snowflake,
            reason=ctx.reason,
        )
        await self.__database_factory.create(text_mute)
        member = ctx.guild.get_member(ctx.member_snowflake)
        if ctx.channel and member:
            try:
                await ctx.channel.set_permissions(
                    target=ctx.member,
                    send_messages=False,
                    add_reactions=False,
                    reason=ctx.reason,
                )
            except discord.Forbidden as e:
                self.__bot.logger.error(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        if not member:
            member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
        self.text_muted_members.update(
            {ctx.member_snowflake: {"name": ctx.display_name}}
        )
        await self.enforce_log(
            author=default_ctx.author,
            channel=ctx.channel,
            duration_value=ctx.duration_value,
            member=member,
            reason=ctx.reason,
            source=source,
        )
        embed = await self.act_embed(ctx=ctx)
        return await state.end(success=embed)

    async def undo_log(self, author, channel, member, source):
        await self.__stream_service.send_log(
            author=author,
            channel=channel,
            identifier="untmute",
            is_modification=True,
            member=member,
            source=source,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="untmute",
            is_modification=True,
            member=member,
        )

    async def undo(self, ctx, default_ctx, source, state):
        await self.__database_factory.delete(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member_snowflake,
        )
        member = ctx.guild.get_member(ctx.member_snowflake)
        if ctx.channel and member:
            try:
                await ctx.channel.set_permissions(
                    target=ctx.member,
                    send_messages=None,
                    add_reactions=None,
                    reason=ctx.reason,
                )
            except discord.Forbidden as e:
                self.__bot.logger.error(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        member = ctx.guild.get_member(ctx.member_snowflake)
        if not member:
            member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
        del self.text_muted_members[ctx.member_snowflake]
        await self.undo_log(
            author=default_ctx.author,
            channel=ctx.channel,
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
            title=f"{self.__emoji.get_random_emoji()} {member_display_name} has been text-muted",
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
            title=f"{self.__emoji.get_random_emoji()} {member_display_name} has been untext-muted",
            description=(f"**User:** {member_str}\n**Channel:** {ctx.channel.mention}"),
            color=discord.Color.yellow(),
        )
        if member:
            embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def migrate(self, kwargs):
        self.__database_factory.update(**kwargs)

    async def is_text_muted(
        self, channel: discord.abc.GuildChannel, member: discord.Member
    ):
        mute = await self.__database_factory.select(
            channel_snowflake=channel.id,
            guild_snowflake=channel.guild.id,
            member_snowflake=member.id,
            singular=True,
        )
        if mute:
            return True
        return False

    async def is_text_muted_then_mute_and_reset_cooldown(
        self, channel: discord.abc.GuildChannel, member: discord.Member
    ):
        if await self.is_text_muted(channel=channel, member=member):
            targets = []
            for target, overwrite in channel.overwrites.items():
                if any(value is not None for value in overwrite._values.values()):
                    if isinstance(target, discord.Member):
                        targets.append(target)
            if member not in targets:
                try:
                    await channel.set_permissions(
                        member,
                        send_messages=False,
                        add_reactions=False,
                        reason="Reinstating active text-mute.",
                    )
                    await self.update_last_text_muted(channel=channel, member=member)
                except discord.Forbidden as e:
                    self.__bot.logger.warning(e)

    async def update_last_text_muted(
        self, channel: discord.abc.GuildChannel, member: discord.Member
    ):
        where_kwargs = {
            "channel_snowflake": channel.id,
            "guild_snowflake": channel.guild.id,
            "member_snowflake": member.id,
        }
        set_kwargs = {
            "last_muted": datetime.now(timezone.utc),
            "reset": False,
        }
        await self.__database_factory.update(
            set_kwargs=set_kwargs,
            where_kwargs=where_kwargs,
        )
        self.__bot.logger.info(
            f"Updated last_muted record for {member.display_name} in {channel.name}."
        )
