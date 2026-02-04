"""text_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the text mute moderation.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.db.alias.alias_service import AliasService
from vyrtuous.db.infractions.tmute.text_mute import TextMute
from vyrtuous.db.mgmt.stream.stream_service import StreamService
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


class TextMuteService(AliasService):

    lines, pages = [], []
    model = TextMute

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        text_mutes = await TextMute.select(singular=False, **where_kwargs)
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
                {
                    "reason": text_mute.reason,
                    "expires_in": DurationObject.from_expires_in(text_mute.expires_in),
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
                TextMuteService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                TextMuteService.pages.extend(
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
        title = f"{get_random_emoji()} Text Mutes {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await TextMuteService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, text_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    TextMuteService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in text_mute_dictionary.get(
                    "text_mutes", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        TextMuteService.lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
                        TextMuteService.lines.append(
                            f"**Expires in:** {channel_dictionary['expires_in']}"
                        )
                        TextMuteService.lines.append(
                            f"**Reason:** {channel_dictionary['reason']}"
                        )
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(TextMuteService.lines),
                            inline=False,
                        )
                        embed = flush_page(
                            embed, TextMuteService.pages, title, guild.name
                        )
                        TextMuteService.lines = []
                        field_count = 0
            if TextMuteService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(TextMuteService.lines),
                    inline=False,
                )
            TextMuteService.pages.append(embed)
        return TextMuteService.pages

    @classmethod
    async def text_mute_overwrite(cls, channel, member):
        kwargs = {
            "channel_snowflake": channel.id,
            "guild_snowflake": channel.guild.id,
            "member_snowflake": member.id,
        }
        text_mute = await TextMute.select(**kwargs, singular=True)
        if text_mute:
            targets = []
            for target, overwrite in channel.overwrites.items():
                if any(value is not None for value in overwrite._values.values()):
                    if isinstance(target, discord.Member):
                        targets.append(target)
            if member not in targets:
                try:
                    await channel.set_permissions(
                        target=member,
                        send_messages=False,
                        add_reactions=False,
                        reason="Reinstating active text-mute overwrite.",
                    )
                    set_kwargs = {
                        "last_muted": datetime.now(timezone.utc),
                        "reset": False,
                    }
                    await TextMute.update(set_kwargs=set_kwargs, where_kwargs=kwargs)
                except discord.Forbidden as e:
                    logger.warning(e)

    @classmethod
    async def enforce(cls, information, message, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        text_mute = TextMute(
            **information["snowflake_kwargs"],
            expires_in=information["expires_in"],
            reason=information["reason"],
        )
        await text_mute.create()
        channel = message.guild.get_channel(
            information["snowflake_kwargs"]["channel_snowflake"]
        )
        if channel:
            try:
                await channel.set_permissions(
                    target=member,
                    send_messages=False,
                    add_reactions=False,
                    reason=information["reason"],
                )
            except discord.Forbidden as e:
                logger.error(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        await StreamService.send_entry(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            duration=information["duration"],
            identifier="tmute",
            member=member,
            message=message,
            reason=information["reason"],
        )
        embed = await TextMuteService.act_embed(information=information)
        return await state.end(success=embed)

    @classmethod
    async def undo(cls, information, message, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        await TextMute.delete(**information["snowflake_kwargs"])
        channel = message.guild.get_channel(
            information["snowflake_kwargs"]["channel_snowflake"]
        )
        if channel:
            try:
                await channel.set_permissions(
                    target=member,
                    send_messages=None,
                    add_reactions=None,
                    reason=information["reason"],
                )
            except discord.Forbidden as e:
                logger.error(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        await StreamService.send_entry(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            identifier="untmute",
            is_modification=True,
            member=member,
            message=message,
        )
        embed = await TextMuteService.undo_embed(information=information)
        return await state.end(success=embed)

    @classmethod
    async def act_embed(cls, information):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Text-Muted",
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
            title=f"{get_random_emoji()} " f"{member.display_name} has been Unmuted",
            description=(
                f"**User:** {member.mention}\n" f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
