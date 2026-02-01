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

    ARGS_MAP = {
        "alias_name": 1,
        "member": 2
    }

    TABLE_NAME = "roles"

    @classmethod
    async def act_embed(cls, information, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        role = guild.get_role(information["snowflake_kwargs"]["role_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name} has been granted a role",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, information, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        role = guild.get_role(information["snowflake_kwargs"]["role_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name}'s role has been revoked",
            description=(
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
        cls, information, message, state
    ):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        added_role = Role(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
            role_snowflake=information["snowflake_kwargs"]["role_snowflake"],
        )
        await added_role.create()

        role = message.guild.get_role(information["snowflake_kwargs"]["role_snowflake"])
        if role:
            await Role.administer_role(
                guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
                member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
                role_snowflake=information["snowflake_kwargs"]["role_snowflake"],
            )

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            duration=information["duration"],
            member=member,
            message=message,
        )

        embed = await RoleAlias.act_embed(
            information=information, source=message
        )

        return await state.end(success=embed)
    

    @classmethod
    async def undo(
        cls, information, message, state
    ):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        await Role.delete(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
            role_snowflake=information["snowflake_kwargs"]["role_snowflake"],
        )

        role = message.guild.get_role(information["snowflake_kwargs"]["role_snowflake"])
        if role:
            await Role.revoke_role(
                guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
                member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
                role_snowflake=information["snowflake_kwargs"]["role_snowflake"],
            )

        await Streaming.send_entry(
            alias=information["alias"],
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            is_modification=True,
            member=member,
            message=message,
        )

        embed = await RoleAlias.undo_embed(
            information=information, source=message
        )

        return await state.end(success=embed)
