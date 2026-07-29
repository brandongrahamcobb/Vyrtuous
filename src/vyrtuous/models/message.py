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
