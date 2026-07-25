"""!/bin/python3
target.py The purpose of this program is to provide the Target properties class.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

ResolvedTarget = Union[
    discord.Guild,
    discord.VoiceChannel,
    discord.StageChannel,
    discord.TextChannel,
    discord.Member,
    discord.Role,
    UUID,
    int,
    str,
]


class TargetObject:

    def __init__(self, target: ResolvedTarget | None = None):
        self.__target = target

    @property
    def target(self) -> ResolvedTarget | None:
        return self.__target

    @target.setter
    def target(self, new_target: ResolvedTarget | None) -> None:
        self.__target = new_target


class Converter(commands.Converter):

    def __init__(self, target_cls=TargetObject):
        self.target_cls = target_cls

    async def convert(
        self, ctx: commands.Context, argument: str
    ) -> TargetObject | None:
        bot: DiscordBot = DiscordBot.get_instance()
        if argument.lower() == "all":
            return self.target_cls("all")
        try:
            return self.target_cls(UUID(argument))
        except ValueError:
            pass
        try:
            id = int(argument)
        except ValueError:
            match = re.search(r"\d+", argument)
            if not match:
                raise commands.BadArgument(
                    f"Could not resolve a valid target ({argument})."
                )
            id = int(match.group())
        guild = bot.get_guild(id)
        if guild:
            return self.target_cls(guild)
        for g in bot.guilds:
            channel = g.get_channel(id)
            if isinstance(
                channel,
                (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
            ):
                return self.target_cls(channel)
            member = g.get_member(id)
            if member:
                return self.target_cls(member)
        guild = ctx.guild
        if isinstance(guild, discord.Guild):
            role = guild.get_role(id)
            if role:
                return self.target_cls(role)
        if id in bot.registry.get(MemberState).active:
            return self.target_cls(id)
        raise commands.BadArgument(f"Could not resolve a valid targeti ({argument}).")


class Transformer(app_commands.Transformer):

    def __init__(self, target_cls=TargetObject):
        self.target_cls = target_cls

    async def transform(
        self, interaction: discord.Interaction, argument: str
    ) -> TargetObject | None:
        bot: DiscordBot = DiscordBot.get_instance()
        if argument.lower() == "all":
            return self.target_cls("all")
        try:
            return self.target_cls(UUID(argument))
        except ValueError:
            pass
        try:
            id = int(argument)
        except ValueError:
            match = re.search(r"\d+", argument)
            if not match:
                raise app_commands.TransformerError(
                    argument, discord.AppCommandOptionType.string, self
                )
            id = int(match.group())
        guild = bot.get_guild(id)
        if guild:
            return self.target_cls(guild)
        for g in bot.guilds:
            channel = g.get_channel(id)
            if isinstance(
                channel,
                (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
            ):
                return self.target_cls(channel)
            member = g.get_member(id)
            if member:
                return self.target_cls(member)
        guild = interaction.guild
        if isinstance(guild, discord.Guild):
            role = guild.get_role(id)
            if role:
                return self.target_cls(role)
        if id in bot.registry.get(MemberState).active:
            return self.target_cls(id)
        raise app_commands.TransformerError(
            argument, discord.AppCommandOptionType.string, self
        )


class Target(Converter):
    def __init__(self):
        super().__init__(TargetObject)


class AppTarget(Transformer):
    def __init__(self):
        super().__init__(TargetObject)
