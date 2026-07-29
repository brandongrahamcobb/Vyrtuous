from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot


class GuildNotFound(commands.GuildNotFound, app_commands.AppCommandError):
    def __init__(self, argument: str):
        commands.GuildNotFound.__init__(self, argument)


class ChannelNotFound(commands.ChannelNotFound, app_commands.AppCommandError):
    def __init__(self, argument: str):
        commands.ChannelNotFound.__init__(self, argument)


class RoleNotFound(commands.RoleNotFound, app_commands.AppCommandError):
    def __init__(self, argument: str):
        commands.RoleNotFound.__init__(self, argument)


class MemberNotFound(commands.MemberNotFound, app_commands.AppCommandError):
    def __init__(self, argument: str):
        commands.MemberNotFound.__init__(self, argument)


class BadArgument(commands.BadArgument, app_commands.AppCommandError):
    def __init__(self, argument: str):
        commands.BadArgument.__init__(self, argument)


class CheckFailure(commands.CheckFailure, app_commands.CheckFailure):
    def __init__(self, message: str):
        commands.CheckFailure.__init__(self, message)


class HasEqualOrLowerRole(app_commands.CheckFailure, commands.CheckFailure):
    def __init__(self, target_rank: str):
        super().__init__(
            message=f"You may not execute this command on this `{target_rank}` because they have an equal or higher role than you in this channel/server."
        )


class ExtensionError(app_commands.AppCommandError, commands.ExtensionError):
    def __init__(self, message: str, name: str):
        commands.ExtensionError.__init__(self, message, name=name)


class NotSysadmin(app_commands.CheckFailure, commands.CheckFailure):
    def __init__(self, message="You lack sufficient permissions of a sysadmin."):
        super().__init__(message)


class NotGuildOwner(commands.CheckFailure):
    def __init__(self, guild_snowflake: int | None):
        name = None
        bot: DiscordBot = DiscordBot.get_instance()
        if guild_snowflake:
            guild = bot.get_guild(guild_snowflake)
            if guild is None:
                raise GuildNotFound(str(guild_snowflake))
            name = guild.name
        message: str = (
            f"You lack sufficient permissions of a guild owner in the requested server{f" ({name})" if name else ""}."
        )
        super().__init__(message)
