"""!/bin/python3
discord_object_service.py The purpose of this program is to provide the DiscordObject module.

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

import asyncio
from typing import Union

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot


class TargetIsBot(commands.CheckFailure):
    def __init__(
        self,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        self._source = ctx or interaction or message
        super().__init__(
            message=f"You cannot execute actions on {self._source.guild.me.mention}."
        )


class DiscordObjectService:
    def to_dict(
        self,
        obj: Union[
            str, discord.abc.GuildChannel, discord.Role, discord.Guild, discord.Member
        ],
    ):
        if not obj:
            return None
        kwargs = {"columns": {}}
        if obj and str(obj).lower() == "all":
            return kwargs
        else:
            kwargs["id"] = obj.id
            kwargs["object"] = obj
            kwargs["type"] = type(obj)
        if isinstance(obj, discord.abc.GuildChannel):
            kwargs.update(
                {
                    "columns": {
                        "channel_snowflake": obj.id,
                        "guild_snowflake": obj.guild.id,
                    },
                    "mention": obj.mention,
                    "name": obj.name,
                }
            )
        elif isinstance(obj, discord.Member):
            kwargs.update(
                {
                    "columns": {
                        "member_snowflake": obj.id,
                        "guild_snowflake": obj.guild.id,
                    },
                    "mention": obj.mention,
                    "name": obj.display_name,
                }
            )
        elif isinstance(obj, discord.Guild):
            kwargs.update(
                {
                    "columns": {"guild_snowflake": obj.id},
                    "id": obj.id,
                    "name": obj.name,
                }
            )
        elif isinstance(obj, discord.Role):
            kwargs.update(
                {
                    "columns": {
                        "guild_snowflake": obj.guild.id,
                        "role_snowflake": obj.id,
                    },
                    "mention": obj.mention,
                    "name": obj.name,
                }
            )
        return kwargs


class MultiConverter(commands.Converter):
    def __init__(self):
        self.__bot = DiscordBot.get_instance()

    async def convert(self, ctx: commands.Context, argument: str):
        if argument and str(argument).lower() == "all":
            return "all"
        channel = commands.VoiceChannelConverter()
        guild = commands.GuildConverter()
        member = commands.MemberConverter()
        role = commands.RoleConverter()
        try:
            channel = await channel.convert(ctx, argument)
            return channel
        except commands.ChannelNotFound as e:
            self.__bot.logger.warning(e)
        try:
            member = await member.convert(ctx, argument)
            return member
        except commands.MemberNotFound as e:
            self.__bot.logger.warning(e)
        try:
            guild = await guild.convert(ctx, argument)
            return guild
        except commands.GuildNotFound as e:
            self.__bot.logger.warning(e)
        try:
            role = await role.convert(ctx, argument)
            return role
        except commands.RoleNotFound as e:
            self.__bot.logger.warning(e)
        raise commands.BadArgument("Argument is not a channel, member, guild, or role.")


class DiscordInteractionObject(app_commands.Transformer):
    async def transform(self, interaction: discord.Interaction, value: str):
        bot = interaction.client
        channel = app_commands.Transform[discord.abc.GuildChannel]
        guild = app_commands.Transform[discord.Guild]
        member = app_commands.Transform[discord.Member]
        role = app_commands.Transform[discord.Role]
        try:
            channel = await channel.transform(interaction, value)
        except app_commands.AppCommandError as e:
            bot.logger.warning(e)
        else:
            return {
                "columns": {
                    "channel_snowflake": channel.id,
                    "guild_snowflake": channel.guild.id,
                },
                "id": channel.id,
                "mention": channel.mention,
                "name": channel.name,
                "type": type(channel),
                "object": channel,
            }
        try:
            member = await member.transform(interaction, value)
        except app_commands.AppCommandError as e:
            bot.logger.warning(e)
        else:
            if interaction.guild and member == interaction.guild.me:
                raise TargetIsBot()
            return {
                "columns": {
                    "guild_snowflake": member.guild.id,
                    "member_snowflake": member.id,
                },
                "id": member.id,
                "mention": member.mention,
                "name": member.display_name,
                "type": type(member),
                "object": member,
            }
        try:
            guild = await guild.transform(interaction, value)
        except app_commands.AppCommandError as e:
            bot.logger.warning(e)
        else:
            return {
                "columns": {"guild_snowflake": guild.id},
                "id": guild.id,
                "name": guild.name,
                "type": type(guild),
                "object": guild,
            }
        try:
            role = await role.transform(interaction, value)
        except app_commands.AppCommandError as e:
            bot.logger.warning(e)
        else:
            return {
                "columns": {
                    "guild_snowflake": role.guild.id,
                    "role_snowflake": role.id,
                },
                "id": role.id,
                "mention": role.mention,
                "name": role.name,
                "type": type(role),
                "object": role,
            }
        raise app_commands.AppCommandError(
            "Argument is not a channel, member, guild, or role."
        )
