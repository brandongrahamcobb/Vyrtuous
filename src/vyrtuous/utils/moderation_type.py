from discord import app_commands
from discord.ext import commands
import discord

class ModerationTypeObject:

    def __init__(self, moderation_type: str):
        self.moderation_type = moderation_type

    @property
    def moderation_type(self) -> str:
        return self._moderation_type

    @moderation_type.setter
    def moderation_type(self, value: str):
        match value:
            case "tmute":
                value = "text_mute"
            case "vmute":
                value = "voice_mute"
            case "untmute":
                value = "untext_mute"
            case "unvmute":
                value = "unvoice_mute"
        self._moderation_type = value

    def __str__(self):
        return self.moderation_type
    
class Converter(commands.Converter):

    def __init__(self, moderation_type_cls = ModerationTypeObject):
        self.moderation_type_cls = moderation_type_cls

    async def convert(self, ctx: commands.Context, arg) -> str:
        return self.moderation_type_cls(arg).moderation_type

class Transformer(app_commands.Transformer):
    
    def __init__(self, moderation_type_cls = ModerationTypeObject):
        self.moderation_type_cls = moderation_type_cls

    async def transform(self, interaction: discord.Interaction, arg) -> str:
        return self.moderation_type_cls(arg).moderation_type

class ModerationType(Converter):
    def __init__(self):
        super().__init__(ModerationTypeObject)

class AppModerationType(Transformer):
    def __init__(self):
        super().__init__(ModerationTypeObject)