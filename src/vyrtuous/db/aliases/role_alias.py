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


import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.roles.role import Role
from vyrtuous.db.mgmt.stream import Streaming
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class RoleAlias(Alias):

    identifier = "role"

    ACT = "role"
    UNDO = "unrole"

    @classmethod
    async def act_embed(cls, infraction_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(infraction_information["infraction_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(infraction_information["infraction_member_snowflake"])
        role = source.guild.get_role(infraction_information["infraction_role_snowflake"])
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
    async def undo_embed(cls, infraction_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(infraction_information["infraction_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(infraction_information["infraction_member_snowflake"])
        role = source.guild.get_role(infraction_information["infraction_role_snowflake"])
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
    async def enforce(
        cls, alias, infraction_information, member, message, state
    ):
        added_role = Role(
            channel_snowflake=infraction_information["infraction_channel_snowflake"],
            guild_snowflake=infraction_information["infraction_guild_snowflake"],
            member_snowflake=infraction_information["infraction_member_snowflake"],
            role_snowflake=infraction_information["infraction_role_snowflake"],
        )
        await added_role.create()

        role = message.guild.get_role(alias.role_snowflake)
        if role:
            await Role.administer_role(
                guild_snowflake=infraction_information["infraction_guild_snowflake"],
                member_snowflake=infraction_information["infraction_member_snowflake"],
                role_snowflake=infraction_information["infraction_role_snowflake"],
            )

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=infraction_information["infraction_channel_snowflake"],
            duration=infraction_information["infraction_duration"],
            is_channel_scope=False,
            is_modification=False,
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await Role.act_embed(
            infraction_information=infraction_information, source=message
        )

        return await state.end(success=embed)
    

    @classmethod
    async def undo(
        cls, alias, infraction_information, member, message, state
    ):
        await Role.delete(
            channel_snowflake=infraction_information["infraction_channel_snowflake"],
            guild_snowflake=infraction_information["infraction_guild_snowflake"],
            member_snowflake=infraction_information["infraction_member_snowflake"],
            role_snowflake=infraction_information["infraction_role_snowflake"],
        )

        role = message.guild.get_role(alias.role_snowflake)
        if role:
            await Role.revoke_role(
                guild_snowflake=infraction_information["infraction_guild_snowflake"],
                member_snowflake=infraction_information["infraction_member_snowflake"],
                role_snowflake=infraction_information["infraction_role_snowflake"],
            )

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=infraction_information["infraction_channel_snowflake"],
            duration="",
            is_channel_scope=False,
            is_modification=infraction_information["infraction_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await Role.undo_embed(
            infraction_information=infraction_information, source=message
        )

        return await state.end(success=embed)
