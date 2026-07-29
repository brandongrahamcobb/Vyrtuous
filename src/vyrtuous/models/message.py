"""!/bin/python3
message.py The purpose of this program is to provide the Message properties class.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.utils.errors.error import BadArgument


class MessageObject:

    def __init__(self, message: discord.Message):
        self.__message = message

    @property
    def message(self) -> discord.Message:
        return self.__message


class Converter(commands.Converter):

    def __init__(self, message_cls=MessageObject):
        self.message_cls = message_cls

    async def convert(self, ctx: commands.Context, argument: str) -> MessageObject:
        try:
            message = await ctx.channel.fetch_message(int(argument))
        except (ValueError, discord.NotFound, discord.HTTPException):
            raise BadArgument("Invalid message ID.")
        return self.message_cls(message)


class Transformer(app_commands.Transformer):

    def __init__(self, message_cls=MessageObject):
        self.message_cls = message_cls

    async def transform(
        self, interaction: discord.Interaction, arg: str
    ) -> MessageObject:
        try:
            if isinstance(
                interaction.channel,
                (discord.VoiceChannel, discord.TextChannel, discord.StageChannel),
            ):
                message = await interaction.channel.fetch_message(int(arg))
            else:
                raise app_commands.TransformerError(arg, self.type, self)
        except (ValueError, discord.NotFound, discord.HTTPException):
            raise app_commands.TransformerError(arg, self.type, self)
        return self.message_cls(message)


class Message(Converter):
    def __init__(self):
        super().__init__(MessageObject)


class AppMessage(Transformer):
    def __init__(self):
        super().__init__(MessageObject)
