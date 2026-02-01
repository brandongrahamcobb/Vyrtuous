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
from vyrtuous.db.aliases.text_mute_alias import TextMuteAlias
from vyrtuous.fields.duration import DurationObject
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_dictionary,
    flush_page,
)
from vyrtuous.utils.logger import logger
from vyrtuous.inc.helpers import CHUNK_SIZE


class TextMute(TextMuteAlias):

    lines, pages = [], []

    PLURAL = "Text Mutes"
    SCOPES = ["channel", "guild", "member"]
    SINGULAR = "Text Mute"
    REQUIRED_ARGS = [
        "channel_snowflake",
        "guild_snowflake",
        "member_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "expires_in",
        "last_muted",
        "reason",
        "reset",
        "updated_at",
    ]

    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
        created_at: datetime = datetime.now(timezone.utc),
        expired: bool = False,
        expires_in: datetime = datetime.now(timezone.utc),
        last_muted: datetime = datetime.now(timezone.utc),
        reason: str = "No reason provided.",
        reset: bool = False,
        updated_at: datetime = datetime.now(timezone.utc),
        **kwargs,
    ):
        self.bot = DiscordBot.get_instance()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.expired = expired
        self.expires_in = expires_in
        self.guild_snowflake = guild_snowflake
        self.last_muted = last_muted
        self.member_snowflake = member_snowflake
        self.member_mention = f"<@{member_snowflake}>" if member_snowflake else None
        self.reason = reason
        self.reset = reset
        self.updated_at = updated_at

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
                TextMute.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                TextMute.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} {TextMute.PLURAL}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await TextMute.build_clean_dictionary(
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
                    TextMute.lines.append(
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
                        TextMute.lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
                        TextMute.lines.append(
                            f"**Expires in:** {channel_dictionary['expires_in']}"
                        )
                        TextMute.lines.append(
                            f"**Reason:** {channel_dictionary['reason']}"
                        )
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(TextMute.lines),
                            inline=False,
                        )
                        embed = flush_page(embed, TextMute.pages, title, guild.name)
                        TextMute.lines = []
                        field_count = 0
            if TextMute.lines:
                embed.add_field(
                    name="Information", value="\n".join(TextMute.lines), inline=False
                )
            TextMute.pages.append(embed)
        return TextMute.pages

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