from copy import copy

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

from datetime import datetime, timedelta, timezone

import discord

from vyrtuous.voice_mute.voice_mute import VoiceMute


class VoiceMuteService:
    __CHUNK_SIZE = 7
    MODEL = VoiceMute

    def __init__(
        self,
        *,
        bot=None,
        database_factory=None,
        dictionary_service=None,
        duration_service=None,
        emoji=None,
        stream_service=None,
    ):
        self.__bot = bot
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL
        self.__dictionary_service = dictionary_service
        self.__duration_service = duration_service
        self.__emoji = emoji
        self.__stream_service = stream_service

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

    async def build_clean_dictionary(self, is_at_home, where_kwargs):
        pages = []
        dictionary = {}
        voice_mutes = await self.__database_factory.select(
            target="user", **where_kwargs
        )
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
                {
                    "reason": voice_mute.reason,
                    "expires_in": self.__duration_service.from_expires_in(
                        voice_mute.expires_in
                    ),
                }
            )
        skipped_guilds = self.__dictionary_service.generate_skipped_guilds(dictionary)
        skipped_members = self.__dictionary_service.generate_skipped_members(dictionary)
        cleaned_dictionary = self.__dictionary_service.clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        # if is_at_home:
        #     if skipped_guilds:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_set_pages(
        #                 skipped=skipped_guilds,
        #                 title="Skipped Servers",
        #             )
        #         )
        #     if skipped_members:
        #         pages.extend(
        #             self.__dictionary_service.generate_skipped_dict_pages(
        #                 skipped=skipped_members,
        #                 title="Skipped Members in Server",
        #             )
        #         )
        return cleaned_dictionary

    async def build_pages(self, object_dict, is_at_home):
        lines, pages = [], []
        title = f"{self.__emoji.get_random_emoji()} Voice Mutes {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await self.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        vmute_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = self.__bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, voice_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in voice_mute_dictionary.get(
                    "voice_mutes", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
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
            pages[0].description = f"**({vmute_n})**"
        return pages

    async def room_mute(self, channel_dict, default_kwargs, reason):
        updated_kwargs = default_kwargs.copy()
        updated_kwargs.update(channel_dict.get("columns", None))
        muted_members, pages, skipped_members, failed_members = [], [], [], []
        guild = self.__bot.get_guild(updated_kwargs.get("guild_snowflake", None))
        for member in channel_dict.get("object", None).members:
            if member.id == updated_kwargs.get("member_snowflake", None):
                continue
            voice_mute = await self.__database_factory.select(
                **updated_kwargs, target="user", singular=True
            )
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_dict.get("id", None):
                    try:
                        await member.edit(mute=True)
                    except Exception as e:
                        self.__bot.logger.warning(
                            f"Unable to voice-mute member "
                            f"{member.display_name} ({member.id}) in channel "
                            f"{channel_dict.get('name', None)} ({channel_dict.get('id', None)}) in guild "
                            f"{guild.name} ({guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
            voice_mute = self.MODEL(
                expires_in=expires_in,
                member_snowflake=member.id,
                reason=reason,
                target="user",
                **updated_kwargs,
            )
            await self.__database_factory.create(voice_mute)
            muted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_dict.get('mention', None)}",
            f"**Muted:** {len(muted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f"**Skipped:** {
                len(channel_dict.get('object', None).members)
                - len(muted_members)
                - len(failed_members)
            }",
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{get_random_emoji()} Room Mute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)
        return pages

    async def room_unmute(self, channel_dict, guild_snowflake):
        unmuted_members, pages, skipped_members, failed_members = [], [], [], []
        guild = self.__bot.get_guild(guild_snowflake)
        where_kwargs = channel_dict.get("columns", None)

        for member in channel_dict.get("object", None).members:
            voice_mute = await self.__database_factory.select(
                target="user", **where_kwargs, member_snowflake=member.id, singular=True
            )
            if not voice_mute:
                skipped_members.append(member)
                continue
            await self.__database_factory.delete(target="user", **where_kwargs)
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_dict.get("id", None):
                    try:
                        await member.edit(mute=False)
                    except Exception as e:
                        self.__bot.logger.warning(
                            f"Unable to undo voice-mute "
                            f"for member {member.display_name} ({member.id}) "
                            f"in channel {channel_dict.get('name', None)} ({channel_dict.get('id', None)}) "
                            f"in guild {guild.name} "
                            f"({guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            unmuted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_dict.get('mention', None)}",
            f"**Unmuted:** {len(unmuted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f"**Skipped:** {
                len(channel_dict.get('object', None).members)
                - len(unmuted_members)
                - len(failed_members)
            }",
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{get_random_emoji()} Room Unmute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)
        return pages

    async def enforce(self, ctx, source, state):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        voice_mute = self.MODEL(
            channel_snowflake=ctx.target_channel_snowflake,
            expires_in=ctx.expires_in,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            reason=ctx.reason,
            target="user",
        )
        await self.__database_factory.create(voice_mute)
        is_channel_scope = False
        if member.voice and member.voice.channel:
            if member.voice.channel.id == ctx.target_channel_snowflake:
                is_channel_scope = True
                try:
                    await member.edit(mute=True, reason=ctx.reason)
                except discord.Forbidden as e:
                    return await state.end(error=str(e).capitalize())
        await self.__stream_service.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            duration=DurationObject.from_expires_in(ctx.expires_in),
            identifier="vmute",
            is_channel_scope=is_channel_scope,
            member=member,
            source=source,
            reason=ctx.reason,
        )
        embed = await self.act_embed(ctx=ctx)
        return await state.end(success=embed)

    async def undo(self, ctx, source, state):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        await self.__database_factory.delete(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
        )
        is_channel_scope = False
        if member.voice and member.voice.channel:
            try:
                is_channel_scope = True
                await member.edit(mute=False)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())
        await self.__stream_service.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="unvmute",
            is_channel_scope=is_channel_scope,
            is_modification=True,
            member=member,
            source=source,
        )
        embed = await self.undo_embed(ctx=ctx)
        return await state.end(success=embed)

    async def act_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        channel = guild.get_channel(ctx.target_channel_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member.display_name} has been voice-muted",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {DurationObject.from_expires_in(ctx.expires_in)}\n"
                f"**Reason:** {ctx.reason}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    async def undo_embed(self, ctx):
        guild = self.__bot.get_guild(ctx.source_guild_snowflake)
        channel = guild.get_channel(ctx.target_channel_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} {member.display_name} has been unmuted",
            description=(
                f"**User:** {member.mention}\n**Channel:** {channel.mention}\n"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
