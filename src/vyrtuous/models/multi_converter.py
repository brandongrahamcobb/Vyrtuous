"""!/bin/python3
discord_object_service.py The purpose of this program is to provide the DiscordObject module.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from typing import Union
from uuid import UUID

import discord
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
        int,
        str,
    ]:
        if argument.lower() == "all":
            return "all"
        bot: DiscordBot = DiscordBot.get_instance()
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
        number = argument.strip("<@>")
        if number.isdigit():
            try:
                member = bot.registry.get(MemberState).active.get(int(number), None)
                if member:
                    return int(number)
            except commands.MemberNotFound as e:
                bot.logger.warning(e)
        raise commands.BadArgument(
            "Argument is not a channel, member, guild, role or UUID."
        )
