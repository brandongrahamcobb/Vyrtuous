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
from vyrtuous.utils.dir_to_classes import skip_db_discovery
from vyrtuous.utils.logger import logger


@skip_db_discovery
class Role(DatabaseFactory):

    ACT = "role"
    CATEGORY = None
    PLURAL = "Roles"
    SCOPES = [None]
    SINGULAR = "Role"
    UNDO = "unrole"

    REQUIRED_INSTANTIATION_ARGS = [
        "guild_snowflake",
        "member_snowflake",
        "role_snowflake",
    ]
    OPTIONAL_ARGS = ["channel_snowflake", "created_at", "updated_at"]

    TABLE_NAME = None

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
