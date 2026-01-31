"""!/bin/python3

aliases.py A discord.py cog containing command aliases for the Vyrtuous bot.

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

from datetime import datetime, timezone

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot

from vyrtuous.db.infractions.ban import Ban
from vyrtuous.db.infractions.flag import Flag
from vyrtuous.db.infractions.role import Role
from vyrtuous.db.infractions.text_mute import TextMute
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.db.infractions.voice_mute import VoiceMute
from vyrtuous.db.mgmt.stream import Streaming
from vyrtuous.utils.invincibility import Invincibility
from vyrtuous.utils.logger import logger


class Aliases(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.alias_help = {
            "ban": [
                "**member**: Tag a member or include their ID",
                "**duration**: m/h/d",
                "**reason**: Reason for ban",
            ],
            "flag": [
                "**member**: Tag a member or include their ID",
                "**reason**: Reason for flag",
            ],
            "role": [
                "**member**: Tag a member or include their ID",
                "**role**: Tag a role or include its ID",
            ],
            "tmute": [
                "**member**: Tag a member or include their ID",
                "**duration**: m/h/d",
                "**reason**: Reason for text-mute",
            ],
            "vegan": ["**member**: Tag a member or include their ID"],
            "vmute": [
                "**member**: Tag a member or include their ID",
                "**duration**: m/h/d",
                "**reason**: Reason for voice-mute",
            ],
        }
        self.category_to_description = {
            "ban": "Toggles a ban.",
            "flag": "Toggles a moderation flag.",
            "role": "Toggles a role to a user.",
            "tmute": "Toggles a mute in text channels.",
            "vegan": "Toggles a vegan.",
            "vmute": "Toggles a mute in voice channels.",
        }
        self.category_to_permission_level = {
            "ban": "Moderator",
            "flag": "Moderator",
            "role": "Coordinator",
            "tmute": "Moderator",
            "vegan": "Moderator",
            "vmute": "Moderator",
        }
        self.bot = bot
        self.invincible_members = Invincibility.get_invincible_members()


    

async def setup(bot: DiscordBot):
    cog = Aliases(bot)
    await bot.add_cog(cog)
