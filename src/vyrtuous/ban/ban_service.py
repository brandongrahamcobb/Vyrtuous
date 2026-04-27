"""!/bin/python3
ban_service.py The purpose of this program is to extend AliasService to service ban infractions.

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

from vyrtuous.ban.ban import Ban


@dataclass
class BanDictionary:
    data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, Any]]]]]] = field(
        default_factory=dict
    )
    skipped_guilds: List[discord.Embed] = field(default_factory=list)
    skipped_members: List[discord.Embed] = field(default_factory=list)


class BanService:
    __CHUNK_SIZE = 12
    MODEL = Ban
    banned_members = {}

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
        self.__bot = bot
        self.__active_member_service = active_member_service
        self.__bot.logger.info(type(self.__active_member_service))
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = Ban
        self.__data_service = data_service
        self.__dictionary_service = dictionary_service
        self.__duration_builder = duration_builder
        self.__emoji = emoji
        self.__stream_service = stream_service

    async def populate(self):
        banned_members = await self.__database_factory.select()
        for banned_member in banned_members:
            guild = self.__bot.get_guild(banned_member.guild_snowflake)
            if not guild:
                continue
            self.banned_members[banned_member.member_snowflake] = {
                "last_active": None,
                "name": banned_member.display_name,
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
            member_snowflake=ctx.member_snowflake,
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
        expired_bans = await self.__database_factory.select(expired=True)
        if expired_bans:
            for expired_ban in expired_bans:
                channel_snowflake = int(expired_ban.channel_snowflake)
                guild_snowflake = int(expired_ban.guild_snowflake)
                member_snowflake = int(expired_ban.member_snowflake)
                kwargs = {
                    "channel_snowflake": channel_snowflake,
                    "guild_snowflake": guild_snowflake,
                    "member_snowflake": member_snowflake,
                }
                guild = self.__bot.get_guild(guild_snowflake)
                if guild is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate guild {guild_snowflake}, cleaning up expired ban."
                    )
                    continue
                channel = guild.get_channel(channel_snowflake)
                if channel is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}, cleaning up expired ban."
                    )
                    continue
                member = guild.get_member(member_snowflake)
                if member is None:
                    await self.__database_factory.delete(**kwargs)
                    self.__bot.logger.info(
                        f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), cleaning up expired ban."
                    )
                    continue
                await self.__database_factory.delete(**kwargs)
                try:
                    await channel.set_permissions(
                        member, view_channel=None, reason="Cleaning up expired ban."
                    )
                except discord.Forbidden as e:
                    self.__bot.logger.error(str(e).capitalize())

    async def toggle_blacklist(self, channel, member_snowflake):
        ban = await self.__database_factory.select(
            channel_snowflake=channel.id,
            member_snowflake=member_snowflake,
            singular=True,
        )
        member = channel.guild.get_member(member_snowflake)
        if member:
            display_name = member.display_name
        else:
            member = self.__active_member_service.active_members.get(
                member_snowflake, None
            )
            if member:
                display_name = member.get("name", None)
            else:
                display_name = "Unknown member"
        if not ban:
            return f"{display_name} is not banned in {channel.mention}."
        where_kwargs = {
            "channel_snowflake": channel.id,
            "member_snowflake": member_snowflake,
        }
        if ban.blacklisted:
            set_kwargs = {"blacklisted": False}
            action = "unlisted"
            await self.__database_factory.update(
                where_kwargs=where_kwargs, set_kwargs=set_kwargs
            )
        else:
            set_kwargs = {"blacklisted": True}
            action = "blacklisted"
            await self.__database_factory.update(
                where_kwargs=where_kwargs, set_kwargs=set_kwargs
            )
        self.__bot.logger.info(
            f"{display_name} ({member_snowflake}) has been ban {action} in {channel.mention}."
        )
        return f"{display_name} ({member_snowflake}) has been ban {action} in {channel.mention}."

    async def clean_overwrites(self):
        bans = await self.__database_factory.select()
        for ban in bans:
            channel_snowflake = int(ban.channel_snowflake)
            guild_snowflake = int(ban.guild_snowflake)
            member_snowflake = int(ban.member_snowflake)
            where_kwargs = {
                "channel_snowflake": channel_snowflake,
                "guild_snowflake": guild_snowflake,
                "member_snowflake": member_snowflake,
            }
            set_kwargs = {"reset": True}
            if (
                not ban.reset
                and ban.last_kicked < datetime.now(timezone.utc) - timedelta(weeks=1)
                and not ban.blacklisted
            ):
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
                        target=member, overwrite=None, reason="Resetting ban overwrite."
                    )
                except discord.Forbidden as e:
                    self.__bot.logger.error(str(e).capitalize())
                await self.__database_factory.update(
                    set_kwargs=set_kwargs, where_kwargs=where_kwargs
                )

    async def build_dictionary(self, obj):
        bans = []
        dictionary = {}
        if isinstance(obj, discord.Guild):
            bans = await self.__database_factory.select(guild_snowflake=obj.id)
        elif isinstance(obj, discord.abc.GuildChannel):
            bans = await self.__database_factory.select(channel_snowflake=obj.id)
        elif isinstance(obj, discord.Member):
            bans = await self.__database_factory.select(member_snowflake=obj.id)
        else:
            bans = await self.__database_factory.select()
        if bans:
            for ban in bans:
                dictionary.setdefault(ban.guild_snowflake, {"members": {}})
                dictionary[ban.guild_snowflake]["members"].setdefault(
                    ban.member_snowflake, {"bans": {}}
                )
                dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                    "bans"
                ].setdefault(ban.channel_snowflake, {})
                dictionary[ban.guild_snowflake]["members"][ban.member_snowflake][
                    "bans"
                ][ban.channel_snowflake] = {
                    "reason": ban.reason,
                    "expires_in": ban.expires_in,
                    "blacklisted": ban.blacklisted,
                }
        return dictionary

    async def build_pages(self, is_at_home, obj):
        lines, pages = [], []
        thumbnail = False

        obj_name = "All Servers"
        if not isinstance(obj, int):
            obj_name = obj.name
        else:
            member = self.__active_member_service.active_member.get(obj, None)
            if member:
                obj_name = member.get("name", None)
            else:
                return "No active bans found."
        title = f"{self.__emoji.get_random_emoji()} Bans for {obj_name}"

        dictionary = await self.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=BanDictionary, dictionary=dictionary
        )

        for guild_snowflake, guild_data in processed_dictionary.data.items():
            ban_n = 0
            field_count = 0
            lines = []
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, ban_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if member:
                    if not isinstance(obj, discord.Member):
                        lines.append(
                            f"**User:** {member.display_name} {member.mention}"
                        )
                    elif not thumbnail:
                        embed.set_thumbnail(url=obj.display_avatar.url)
                        thumbnail = True
                else:
                    display_name = self.__active_member_service.active_members.get(
                        member_snowflake, None
                    )
                    if member:
                        display_name = member.get("name", None)
                        lines.append(f"**User:** {display_name} ({member_snowflake})")
                for channel_snowflake, channel_dictionary in ban_dictionary.get(
                    "bans"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(obj, discord.abc.GuildChannel):
                        lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(obj, discord.Member):
                        lines.append(
                            f"**Expires in:** {self.__duration_builder.from_timestamp(channel_dictionary['expires_in']).to_unix_ts()}"
                        )
                        lines.append(f"**Reason:** {channel_dictionary['reason']}")
                        lines.append(
                            f"**Blacklisted:** {channel_dictionary['blacklisted']}"
                        )
                    ban_n += 1
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
                    name="Information", value="\n".join(lines), inline=False
                )
            original_description = embed.description or ""
            embed.description = f"**{original_description} ({ban_n})**"
            pages.append(embed)
        if is_at_home:
            pages.extend(processed_dictionary.skipped_guilds)
            pages.extend(processed_dictionary.skipped_members)
        return pages

    async def build_blacklist_pages(self, is_at_home, obj):
        lines, pages = [], []
        thumbnail = False

        obj_name = "All Servers"
        if not isinstance(obj, int):
            obj_name = obj.name
        else:
            member = self.__active_member_service.active_member.get(obj, None)
            if member:
                obj_name = member.get("name", None)
            else:
                return "No active bans found."
        title = f"{self.__emoji.get_random_emoji()} Blacklists for {obj_name}"

        dictionary = await self.build_dictionary(obj=obj)
        processed_dictionary = await self.__dictionary_service.process_dictionary(
            cls=BanDictionary, dictionary=dictionary
        )

        for guild_snowflake, guild_data in processed_dictionary.data.items():
            ban_n = 0
            field_count = 0
            lines = []
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, ban_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if member:
                    if not isinstance(obj, discord.Member):
                        lines.append(
                            f"**User:** {member.display_name} {member.mention}"
                        )
                    elif not thumbnail:
                        embed.set_thumbnail(url=obj.display_avatar.url)
                        thumbnail = True
                else:
                    display_name = self.__active_member_service.active_members.get(
                        member_snowflake, None
                    )
                    if member:
                        display_name = member.get("name", None)
                        lines.append(f"**User:** {display_name} ({member_snowflake})")
                for channel_snowflake, channel_dictionary in ban_dictionary.get(
                    "bans"
                ).items():
                    if not channel_dictionary["blacklisted"]:
                        continue
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(obj, discord.abc.GuildChannel):
                        lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(obj, discord.Member):
                        lines.append(
                            f"**Blacklisted:** {channel_dictionary['blacklisted']}"
                        )
                    ban_n += 1
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
                    name="Information", value="\n".join(lines), inline=False
                )
            original_description = embed.description or ""
            embed.description = f"**{original_description} ({ban_n})**"
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
            if not guild:
                continue
            channel = guild.get_channel(obj.channel_snowflake)
            if not channel:
                continue
            member = guild.get_member(obj.member_snowflake)
            if channel and member:
                try:
                    await channel.set_permissions(target=member, view_channel=None)
                except discord.Forbidden as e:
                    self.__bot.logger.error(str(e).capitalize())
            if not member:
                member = self.__active_member_service.active_members.get(
                    obj.member_snowflake, None
                )
            await self.undo_log(
                author=author,
                channel=channel,
                member=member,
                source=source,
            )

    async def enforce_log(
        self, author, channel, duration_value, is_channel_scope, member, reason, source
    ):
        await self.__stream_service.send_log(
            author=author,
            channel=channel,
            duration_value=duration_value,
            identifier="ban",
            is_channel_scope=is_channel_scope,
            member=member,
            reason=reason,
            source=source,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="ban",
            duration_value=duration_value,
            reason=reason,
            member=member,
        )

    async def enforce(self, ctx, default_ctx, source, state):
        ban = self.MODEL(
            channel_snowflake=ctx.channel.id,
            display_name=ctx.display_name,
            expires_in=ctx.expires_in,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member_snowflake,
            reason=ctx.reason,
        )
        await self.__database_factory.create(ban)
        is_channel_scope = False
        member = ctx.guild.get_member(ctx.member_snowflake)
        if ctx.channel and member:
            try:
                await ctx.channel.set_permissions(
                    member,
                    view_channel=False,
                    reason=ctx.reason,
                )
                if (
                    member.voice
                    and member.voice.channel
                    and member.voice.channel.id == ctx.channel.id
                ):
                    is_channel_scope = True
                    await member.move_to(None, reason=ctx.reason)
                    where_kwargs = {
                        "channel_snowflake": ctx.channel.id,
                        "guild_snowflake": ctx.guild.id,
                        "member_snowflake": ctx.member_snowflake,
                    }
                    set_kwargs = {"last_kicked": datetime.now(timezone.utc)}
                    await self.__database_factory.update(
                        set_kwargs=set_kwargs,
                        where_kwargs=where_kwargs,
                    )
            except discord.Forbidden as e:
                self.__bot.logger.error(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        if not member:
            member = self.__active_member_service.active_members.get(
                ctx.member_snowflake, None
            )
        await self.enforce_log(
            author=default_ctx.author,
            channel=ctx.channel,
            duration_value=ctx.duration_value,
            is_channel_scope=is_channel_scope,
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
            identifier="unban",
            is_modification=True,
            member=member,
            source=source,
        )
        await self.__data_service.save_data(
            author=author,
            channel=channel,
            identifier="unban",
            is_modification=True,
            member=member,
        )

    async def undo(self, ctx, default_ctx, source, state):
        ban = await self.__database_factory.select(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.member_snowflake,
            singular=True,
        )
        if not ban.blacklisted:
            await self.__database_factory.delete(
                channel_snowflake=ctx.channel.id,
                guild_snowflake=ctx.guild.id,
                member_snowflake=ctx.member_snowflake,
            )
            member = ctx.guild.get_member(ctx.member_snowflake)
            if ctx.channel and member:
                try:
                    await ctx.channel.set_permissions(member, view_channel=None)
                except discord.Forbidden as e:
                    self.__bot.logger.error(str(e).capitalize())
                    return await state.end(error=str(e).capitalize())
            if not member:
                member = self.__active_member_service.active_members.get(
                    ctx.member_snowflake, None
                )
            await self.undo_log(
                author=default_ctx.author,
                channel=ctx.channel,
                member=member,
                source=source,
            )
            embed = await self.undo_embed(ctx=ctx)
        else:
            embed = await self.blacklisted_block_embed(ctx=ctx)
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
            title=f"{self.__emoji.get_random_emoji()} {member_display_name} has been banned",
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
            title=f"{self.__emoji.get_random_emoji()} {member_display_name} has been unbanned",
            description=(f"**User:** {member_str}\n**Channel:** {ctx.channel.mention}"),
            color=discord.Color.yellow(),
        )
        if member:
            embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def blacklisted_block_embed(self, ctx):
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
            title=f"{self.__emoji.get_random_emoji()} {member_display_name} is blacklisted",
            description=f"**User:** {member_str}\n**Channel:** {ctx.channel.mention}\nUse {self.__bot.config['discord_command_prefix']}blacklist to unblock.",
            color=discord.Color.yellow(),
        )
        if member:
            embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def migrate(self, kwargs):
        await self.__database_factory.update(**kwargs)

    async def is_banned(
        self, channel: discord.abc.GuildChannel, member: discord.Member
    ):
        ban = await self.__database_factory.select(
            channel_snowflake=channel.id, member_snowflake=member.id
        )
        if ban:
            return True
        return False

    async def is_banned_then_kick_and_reset_cooldown(
        self, channel: discord.abc.GuildChannel, member: discord.Member
    ):
        if await self.is_banned(channel=channel, member=member):
            if (
                member.voice
                and member.voice.channel
                and member.voice.channel.id == channel.id
            ):
                await self.kick(channel=channel, member=member)
            targets = []
            for target, overwrite in channel.overwrites.items():
                if any(value is not None for value in overwrite._values.values()):
                    if isinstance(target, discord.Member):
                        targets.append(target)
            if member not in targets:
                await self.toggle_view_channel(
                    channel=channel, member=member, view_channel=False
                )

    async def toggle_view_channel(
        self,
        channel: discord.abc.GuildChannel,
        member: discord.Member,
        view_channel: bool,
    ):
        try:
            await channel.set_permissions(
                member,
                view_channel=view_channel,
                reason=f"Toggled ban {'off' if not view_channel else 'on'}.",
            )
        except discord.Forbidden as e:
            self.__bot.logger.warning(e)

    async def kick(self, channel: discord.abc.GuildChannel, member: discord.Member):
        try:
            await member.move_to(None, reason="Reinstating active ban.")
            await self.update_last_kicked(channel=channel, member=member)
        except discord.Forbidden as e:
            self.__bot.logger.warning(e)

    async def update_last_kicked(
        self, channel: discord.abc.GuildChannel, member: discord.Member
    ):
        where_kwargs = {"channel_snowflake": channel.id, "member_snowflake": member.id}
        set_kwargs = {
            "last_kicked": datetime.now(timezone.utc),
            "reset": False,
        }
        await self.__database_factory.update(
            set_kwargs=set_kwargs,
            where_kwargs=where_kwargs,
        )
        self.__bot.logger.info(
            f"Updated last_kicked record for {member.display_name} in {channel.name}."
        )
