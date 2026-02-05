import discord
from discord.ext import commands


class InsufficientPermissions(commands.CheckFailure):
    def __init__(self, message="Member has insufficient permissions."):
        super().__init__(message)


class HasEqualOrLowerRole(InsufficientPermissions):
    def __init__(self, target_rank=str):
        super().__init__(
            message=f"You may not execute this command on this `{target_rank}` because they have equal or higher role than you in this channel/server."
        )


class NotAdministrator(InsufficientPermissions):
    def __init__(
        self,
        message="Member is not an administrator in this server.",
    ):
        super().__init__(message)


class NotCoordinator(InsufficientPermissions):
    def __init__(
        self,
        message="Member is not a coordinator in this channel.",
    ):
        super().__init__(message)


class NotDeveloper(InsufficientPermissions):
    def __init__(self, message="Member is not a developer."):
        super().__init__(message)


class NotModerator(InsufficientPermissions):
    def __init__(
        self,
        message="Member is not a moderator in this channel.",
    ):
        super().__init__(message)


class NotGuildOwner(InsufficientPermissions):
    def __init__(
        self,
        message="Member is not a guild owner in this server.",
    ):
        super().__init__(message)


class NotSysadmin(InsufficientPermissions):
    def __init__(self, message="Member is not a sysadmin."):
        super().__init__(message)


class DiscordObjectNotFound(commands.CheckFailure):
    "Returns an error if a channel, guild, member or role is not found."

    def __init__(self, target: str, *, message: str | None = None):
        super().__init__(
            message=message
            or f"Unable to resolve a valid channel, guild, member or role for target (`{target}`)."
        )


class GuildChannelNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            message=f"Unable to resolve a valid Discord guild channel with the provided context and target (`{target}`).",
            target=target,
        )


class GuildNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            message=f"Unable to resolve a valid Discord guild with the provided context and target (`{target}`).",
            target=target,
        )


class GuildMemberNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            message=f"Unable to resolve a valid Discord guild member with the provided context and target (`{target}`).",
            target=target,
        )


class GuildRoleNotFound(DiscordObjectNotFound):

    def __init__(self, target: str):
        super().__init__(
            message=f"Unable to resolve a valid Discord guild role with the provided context and target (`{target}`).",
            target=target,
        )


class TargetIsBot(commands.CheckFailure):
    def __init__(
        self,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        if (ctx is None) == (interaction is None) == (message is None):
            raise DiscordSourceNotFound()
        self._source = ctx or interaction or message
        super().__init__(
            message=f"You cannot execute actions on {self._source.guild.me.mention}."
        )


class NotAtHome(commands.CheckFailure):
    def __init__(self, message="You are not in the home server and cannot do this."):
        super().__init__(message)


class DiscordSourceNotFound(commands.CheckFailure):

    def __init__(self):
        super().__init__(
            message="Unable to resolve a valid Discord object due to missing ctx (commands.Context), interaction (discord.Interaction) or message (discord.Message)."
        )
