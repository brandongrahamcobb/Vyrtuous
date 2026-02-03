"""voice_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the voice mute moderation.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.base.alias.alias_service import AliasService
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute
from vyrtuous.db.mgmt.stream.stream_service import StreamService
from vyrtuous.fields.duration import DurationObject
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class VoiceMuteService(AliasService):

    lines, pages = [], []
    model = VoiceMute

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        voice_mutes = await VoiceMute.select(target="user", **where_kwargs)
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
                    "expires_in": DurationObject.from_expires_in(voice_mute.expires_in),
                }
            )
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                VoiceMuteService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                VoiceMuteService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Voice Mutes {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await VoiceMuteService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, voice_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    VoiceMuteService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
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
                        VoiceMuteService.lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
                        VoiceMuteService.lines.append(
                            f"**Expires in:** {channel_dictionary['expires_in']}"
                        )
                        VoiceMuteService.lines.append(
                            f"**Reason:** {channel_dictionary['reason']}"
                        )
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(VoiceMuteService.lines),
                            inline=False,
                        )
                        embed = flush_page(
                            embed, VoiceMuteService.pages, title, guild.name
                        )
                        VoiceMuteService.lines = []
                        field_count = 0
            if VoiceMuteService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(VoiceMuteService.lines),
                    inline=False,
                )
            VoiceMuteService.pages.append(embed)
        return VoiceMuteService.pages

    @classmethod
    async def room_mute(cls, channel_dict, guild_snowflake, reason, snowflake_kwargs):
        bot = DiscordBot.get_instance()
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        member_snowflake = snowflake_kwargs.get("member_snowflake", None)
        muted_members, pages, skipped_members, failed_members = [], [], [], []
        guild = bot.get_guild(guild_snowflake)
        where_kwargs = channel_dict.get("columns", None)

        for member in channel_dict.get("object", None).members:
            if member.id == member_snowflake:
                continue
            voice_mute = await VoiceMute.select(
                **where_kwargs, target="user", singular=True
            )
            if voice_mute:
                skipped_members.append(member)
                continue
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_dict.get("id", None):
                    try:
                        await member.edit(mute=True)
                    except Exception as e:
                        logger.warning(
                            f"Unable to voice-mute member "
                            f"{member.display_name} ({member.id}) in channel "
                            f"{channel_dict.get('name', None)} ({channel_dict.get('id', None)}) in guild "
                            f"{guild.name} ({guild.id}). "
                            f"{str(e).capitalize()}"
                        )
                        failed_members.append(member)
            expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
            voice_mute = VoiceMute(
                expires_in=expires_in,
                member_snowflake=member.id,
                reason=reason,
                target="user",
                **where_kwargs,
            )
            await voice_mute.create()
            muted_members.append(member)
        description_lines = [
            f"**Channel:** {channel_dict.get("mention", None)}",
            f"**Muted:** {len(muted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_dict.get('object', None).members) \
                - len(muted_members) \
                - len(failed_members)
            }',
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{get_random_emoji()} Room Mute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)
        return pages

    @classmethod
    async def room_unmute(cls, channel_dict, guild_snowflake):
        unmuted_members, pages, skipped_members, failed_members = [], [], [], []
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        where_kwargs = channel_dict.get("columns", None)

        for member in channel_dict.get("object", None).members:
            voice_mute = await VoiceMute.select(
                target="user", **where_kwargs, member_snowflake=member.id, singular=True
            )
            if not voice_mute:
                skipped_members.append(member)
                continue
            await VoiceMute.delete(target="user", **where_kwargs)
            if member.voice and member.voice.channel:
                if member.voice.channel.id == channel_dict.get("id", None):
                    try:
                        await member.edit(mute=False)
                    except Exception as e:
                        logger.warning(
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
            f"**Channel:** {channel_dict.get("mention", None)}",
            f"**Unmuted:** {len(unmuted_members)} users",
            f"**Failed:** {len(failed_members)} users",
            f'**Skipped:** {len(channel_dict.get('object', None).members) \
                - len(unmuted_members) \
                - len(failed_members)
            }',
        ]
        embed = discord.Embed(
            description="\n".join(description_lines),
            title=f"{get_random_emoji()} Room Unmute Summary",
            color=discord.Color.blurple(),
        )
        pages.append(embed)
        return pages

    @classmethod
    async def enforce(cls, information, message, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        voice_mute = VoiceMute(
            **information["snowflake_kwargs"],
            expires_in=information["expires_in"],
            reason=information["reason"],
            target="user",
        )
        await voice_mute.create()
        is_channel_scope = False
        if member.voice and member.voice.channel:
            if (
                member.voice.channel.id
                == information["snowflake_kwargs"]["channel_snowflake"]
            ):
                is_channel_scope = True
                try:
                    await member.edit(mute=True, reason=information["reason"])
                except discord.Forbidden as e:
                    return await state.end(error=str(e).capitalize())
        await StreamService.send_entry(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            duration=information["duration"],
            identifier="vmute",
            is_channel_scope=is_channel_scope,
            member=member,
            message=message,
            reason=information["reason"],
        )
        embed = await VoiceMuteService.act_embed(information=information)
        return await state.end(success=embed)

    @classmethod
    async def undo(cls, information, message, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        await VoiceMute.delete(**information["snowflake_kwargs"])
        is_channel_scope = False
        if member.voice and member.voice.channel:
            try:
                is_channel_scope = True
                await member.edit(mute=False)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())
        await StreamService.send_entry(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            identifier="unvmute",
            is_channel_scope=is_channel_scope,
            is_modification=True,
            member=member,
            message=message,
        )
        embed = await VoiceMuteService.undo_embed(information=information)
        return await state.end(success=embed)

    @classmethod
    async def act_embed(cls, information):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name} has been voice-muted",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {information['duration']}\n"
                f"**Reason:** {information['reason']}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, information):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been unmuted",
            description=(
                f"**User:** {member.mention}\n" f"**Channel:** {channel.mention}\n"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
