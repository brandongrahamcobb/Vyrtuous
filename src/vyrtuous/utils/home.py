"""home.py The purpose of this program is to provide the at_home utility module.

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
