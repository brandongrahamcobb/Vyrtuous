
''' discord_message_service.py  The purpose of this program is to handle messages in Discord.

    Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
'''
from discord.ext import commands
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.snowflake import *
import asyncio
import discord

class Paginator(discord.ui.View):
    def __init__(self, bot, interaction, pages, *, info_view=None, timeout=60):
        super().__init__(timeout=timeout)
        self.bot = bot
        self.interaction = interaction
        self.pages = pages
        self.current_page = 0
        self.info_view = info_view

    async def start(self):
        view = self.combine_views(self.info_view) if self.info_view else self
        if isinstance(self.interaction, discord.Interaction):
            await self.interaction.response.defer()
            self.message = await self.interaction.followup.send(
                embed=self.pages[self.current_page],
                view=view
            )
        else:  # commands.Context
            self.message = await self.interaction.send(
                embed=self.pages[self.current_page],
                view=view
            )
        return self.message

    @discord.ui.button(label="Previous", style=discord.ButtonStyle.primary)
    async def prev_button(self, interaction: discord.Interaction, button: discord.ui.Button):
        self.current_page = max(self.current_page - 1, 0)
        await interaction.response.edit_message(embed=self.pages[self.current_page], view=self.combine_views(self.info_view))

    @discord.ui.button(label="Next", style=discord.ButtonStyle.primary)
    async def next_button(self, interaction: discord.Interaction, button: discord.ui.Button):
        self.current_page = min(self.current_page + 1, len(self.pages) - 1)
        await interaction.response.edit_message(embed=self.pages[self.current_page], view=self.combine_views(self.info_view))

    def combine_views(self, extra_view: discord.ui.View = None):
        combined = discord.ui.View(timeout=self.timeout)
        combined.add_item(self.prev_button)
        combined.add_item(self.next_button)
        if extra_view:
            for item in extra_view.children:
                combined.add_item(item)  # info/state buttons added after navigation
        return combined
