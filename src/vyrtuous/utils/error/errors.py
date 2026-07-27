from discord import app_commands
from discord.ext import commands


class GuildNotFound(GuildNotFound, app_commands.AppCommandError):
    def __init__(self, argument: str):
        GuildNotFound.__init__(self, argument)


class ChannelNotFound(commands.ChannelNotFound, app_commands.AppCommandError):
    def __init__(self, argument: str):
        commands.ChannelNotFound.__init__(self, argument)


class CheckFailure(commands.CheckFailure, app_commands.CheckFailure):
    def __init__(self, message: str):
        commands.CheckFailure.__init__(self, message)
