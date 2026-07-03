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

from uuid import UUID

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState


class TargetIsBot(commands.CheckFailure):
    def __init__(
        self,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        self._source = ctx or interaction or message
        if self._source is None or self._source.guild is None:
            return
        super().__init__(
            message=f"You cannot execute actions on {self._source.guild.me.mention}."
        )


class MultiConverter(commands.Converter):
    def __init__(self):
        self.__bot = DiscordBot.get_instance()

    async def convert(self, ctx: commands.Context, argument: str):
        if argument and str(argument).lower() == "all":
            return None
        try:
            uuid = UUID(argument)
            return uuid
        except ValueError as e:
            self.__bot.logger.warning(e)
        voice_channel = commands.VoiceChannelConverter()
        stage_channel = commands.StageChannelConverter()
        guild = commands.GuildConverter()
        member = commands.MemberConverter()
        role = commands.RoleConverter()
        try:
            voice_channel = await voice_channel.convert(ctx, argument)
            return voice_channel
        except commands.ChannelNotFound as e:
            self.__bot.logger.warning(e)
        try:
            stage_channel = await stage_channel.convert(ctx, argument)
            return stage_channel
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
        if isinstance(argument, int):
            try:
                member = self.__bot.registry.get(MemberState).active.get(argument, None)
                if member:
                    return argument
            except commands.MemberNotFound as e:
                self.__bot.logger.warning(e)
        raise commands.BadArgument(
            "Argument is not a channel, member, guild, role or UUID."
        )
