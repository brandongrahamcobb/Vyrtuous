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

import re
from typing import Union
from uuid import UUID

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState


class MultiConverter(commands.Converter):

    async def convert(self, ctx: commands.Context, argument: str) -> Union[
        discord.VoiceChannel,
        discord.StageChannel,
        discord.Guild,
        discord.Member,
        discord.Role,
        UUID,
        str,
        None,
    ]:
        bot: DiscordBot = DiscordBot.get_instance()
        if argument and str(argument).lower() == "all":
            return argument
        try:
            uuid = UUID(argument)
            return uuid
        except ValueError as e:
            bot.logger.warning(e)
        voice_channel_converter = commands.VoiceChannelConverter()
        stage_channel_converter = commands.StageChannelConverter()
        guild_converter = commands.GuildConverter()
        member_converter = commands.MemberConverter()
        role_converter = commands.RoleConverter()
        try:
            voice_channel = await voice_channel_converter.convert(ctx, argument)
            return voice_channel
        except commands.ChannelNotFound as e:
            bot.logger.warning(e)
        try:
            stage_channel = await stage_channel_converter.convert(ctx, argument)
            return stage_channel
        except commands.ChannelNotFound as e:
            bot.logger.warning(e)
        try:
            member = await member_converter.convert(ctx, argument)
            return member
        except commands.MemberNotFound as e:
            bot.logger.warning(e)
        try:
            guild = await guild_converter.convert(ctx, argument)
            return guild
        except commands.GuildNotFound as e:
            bot.logger.warning(e)
        try:
            role = await role_converter.convert(ctx, argument)
            return role
        except commands.RoleNotFound as e:
            bot.logger.warning(e)
        if isinstance(argument, int):
            try:
                member = bot.registry.get(MemberState).active.get(argument, None)
                if member:
                    return argument
            except commands.MemberNotFound as e:
                bot.logger.warning(e)
        raise commands.BadArgument(
            "Argument is not a channel, member, guild, role or UUID."
        )


def transform(interaction: discord.Interaction, argument: str | None) -> Union[
    discord.VoiceChannel,
    discord.StageChannel,
    discord.TextChannel,
    discord.Guild,
    discord.Member,
    discord.Role,
    UUID,
    str,
    int,
    None,
]:
    id = None
    bot: DiscordBot = DiscordBot.get_instance()
    if argument is None:
        return None
    if argument and str(argument).lower() == "all":
        return argument
    try:
        uuid = UUID(argument)
        return uuid
    except ValueError as e:
        bot.logger.warning(e)
    try:
        id = int(argument)
    except ValueError:
        try:
            match = re.search(r"\d+", argument)
            if match:
                id = int(match.group())
        except ValueError:
            raise app_commands.CheckFailure("Value is not a valid integer or mention.")
    if id:
        channel = bot.get_channel(id)
        if isinstance(
            channel, (discord.VoiceChannel, discord.StageChannel, discord.TextChannel)
        ):
            return channel
        guild = bot.get_guild(id)
        if guild:
            return guild
        else:
            guild = interaction.guild
            if guild is None:
                raise app_commands.CheckFailure(
                    "This command must be executed in a server."
                )
        member = guild.get_member(id)
        if member:
            return member
        role = guild.get_role(id)
        if role:
            return role
        try:
            member = bot.registry.get(MemberState).active.get(id, None)
            if member:
                return id
        except commands.MemberNotFound as e:
            bot.logger.warning(e)
    raise commands.BadArgument(
        "Argument is not a channel, member, guild, role or UUID."
    )
