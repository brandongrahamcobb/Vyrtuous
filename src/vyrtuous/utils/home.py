from typing import Union

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot


class NotAtHome(commands.CheckFailure):
    def __init__(self, message="You are not in the home server and cannot do this."):
        super().__init__(message)


def at_home(
    source: Union[commands.Context, discord.Interaction, discord.Message],
) -> bool:
    bot = DiscordBot.get_instance()
    if source.guild.id == int(bot.config["discord_testing_guild_snowflake"]):
        return True
    return False
