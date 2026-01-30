from discord.ext import commands


class DiscordSourceNotFound(commands.CheckFailure):

    def __init__(self):
        super().__init__(
            message="Unable to resolve a valid Discord object due to missing ctx (commands.Context), interaction (discord.Interaction) or message (discord.Message)."
        )
