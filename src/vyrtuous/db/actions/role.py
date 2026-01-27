"""flag.py The purpose of this program is to inherit from DatabaseFactory to provide the flag moderation.

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

from datetime import datetime
from typing import Optional

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    clean_guild_dictionary,
    flush_page,
)


class Role(DatabaseFactory):

    ACT = "role"
    CATEGORY = None
    PLURAL = "Roles"
    SCOPES = ["channel", "member"]
    SINGULAR = "Role"
    UNDO = "unrole"

    REQUIRED_INSTANTIATION_ARGS = [
        "guild_snowflake",
        "member_snowflake",
        "role_snowflake",
    ]
    OPTIONAL_ARGS = ["channel_snowflake", "created_at", "updated_at"]

    TABLE_NAME = "roles"

    def __init__(
        self,
        guild_snowflake: int,
        member_snowflake: int,
        role_snowflake: int,
        channel_snowflake: Optional[int] = None,
        created_at: Optional[datetime] = None,
        updated_at: Optional[datetime] = None,
        **kwargs,
    ):
        super().__init__()
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>" if channel_snowflake else None
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.member_snowflake = member_snowflake
        self.role_snowflake = role_snowflake
        self.updated_at = updated_at

    @classmethod
    async def act_embed(cls, action_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(action_information["action_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        role = source.guild.get_member(action_information["action_role_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name} has been granted a role",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, action_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(action_information["action_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(action_information["action_member_snowflake"])
        role = source.guild.get_member(action_information["action_role_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name}'s role has been revoked",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def administer_role(cls, guild_snowflake, member_snowflake, role_snowflake):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        role = guild.get_role(role_snowflake)
        try:
            await member.add_roles(role, reason="Granting role.")
        except discord.Forbidden as e:
            logger.error(str(e).capitalize())

    @classmethod
    async def revoke_role(cls, guild_snowflake, member_snowflake, role_snowflake):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        role = guild.get_role(role_snowflake)
        try:
            await member.remove_roles(role, reason="Revoking role.")
        except discord.Forbidden as e:
            logger.error(str(e).capitalize())

    @classmethod
    async def added_role(
        cls,
        category_class,
        category_role_class,
        guild_snowflake,
        member_snowflake,
        role_snowflake,
    ):
        kwargs = {"guild_snowflake": guild_snowflake, "role_snowflake": role_snowflake}
        role = await category_role_class.select(singular=True, **kwargs)
        if role:
            if hasattr(role, "channel_snowflake"):
                kwargs.update({"channel_snowflake": role.channel_snowflake})
                msg = f"Member ({member_snowflake}) was granted the role ({role_snowflake}) for category ({category_class.__name__()}) related to channel ({role.channel_snowflake}) in guild ({guild_snowflake})."
            else:
                msg = f"Member ({member_snowflake}) was granted the role ({role_snowflake}) for category ({category_class.__name__()}) in guild ({guild_snowflake})."
            category = category_class(**kwargs)
            await category.create()
            logger.info(msg)
        else:
            return

    @classmethod
    async def removed_role(
        cls,
        category_class,
        category_role_class,
        guild_snowflake,
        member_snowflake,
        role_snowflake,
    ):
        kwargs = {"guild_snowflake": guild_snowflake, "role_snowflake": role_snowflake}
        role = await category_role_class.select(singular=True, **kwargs)
        if role:
            if hasattr(role, "channel_snowflake"):
                kwargs.update({"channel_snowflake": role.channel_snowflake})
                msg = f"Member ({member_snowflake}) was revoked the role ({role_snowflake}) for category ({category_class.__name__()}) related to channel ({role.channel_snowflake}) in guild ({guild_snowflake})."
            else:
                msg = f"Member ({member_snowflake}) was revoked the role ({role_snowflake}) for category ({category_class.__name__()}) in guild ({guild_snowflake})."
            await category_class.delete(**kwargs)
            logger.info(msg)
        else:
            return

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        thumbnail = False
        kwargs = object_dict.get("columns", None)
        title = f"{get_random_emoji()} {Role.PLURAL} {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        roles = await Role.select(**kwargs)

        for role in roles:
            guild_dictionary.setdefault(role.guild_snowflake, {"members": {}})
            guild_dictionary[role.guild_snowflake]["members"].setdefault(
                role.member_snowflake, {"roles": {}}
            )
            guild_dictionary[role.guild_snowflake]["members"][role.member_snowflake][
                "roles"
            ].setdefault(role.channel_snowflake, {})
            guild_dictionary[role.guild_snowflake]["members"][role.member_snowflake][
                "roles"
            ][role.channel_snowflake] = {"role": role.role_snowflake}

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_members = generate_skipped_members(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, role_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    lines.append(f"**User:** {member.display_name} {member.mention}")
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in role_dictionary.get(
                    "roles"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    role = channel_dictionary.get(role, None)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        lines.append(f"**Channel:** {channel.mention}")
                    lines.append(f"**Role:** {role.mention}")
                    field_count += 1
                    if field_count >= chunk_size:
                        embed.add_field(
                            name="Information", value="\n".join(lines), inline=False
                        )
                        embed, field_count = flush_page(embed, pages, title, guild.name)
                        lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_members:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_members,
                    title="Skipped Members in Server",
                )
        return pages
