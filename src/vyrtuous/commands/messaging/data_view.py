"""moderation_view.py The purpose of this program is to provide the moderation view utility class.

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
from io import BytesIO
from pathlib import Path

import discord
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from vyrtuous.base.infraction import Infraction
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.db.infractions.ban.ban import Ban
from vyrtuous.db.infractions.flag.flag import Flag
from vyrtuous.db.infractions.tmute.text_mute import TextMute
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute
from vyrtuous.db.roles.admin.administrator import Administrator
from vyrtuous.db.roles.coord.coordinator import Coordinator
from vyrtuous.db.roles.dev.developer import Developer
from vyrtuous.db.roles.mod.moderator import Moderator
from vyrtuous.db.roles.owner.guild_owner_service import (
    NotGuildOwner,
    is_guild_owner_wrapper,
)
from vyrtuous.db.roles.sysadmin.sysadmin_service import (
    NotSysadmin,
    is_sysadmin_wrapper,
)
from vyrtuous.utils.dir_to_classes import dir_to_classes


class DataView(discord.ui.View):

    def __init__(
        self, interaction: discord.Interaction
    ):
        super().__init__(timeout=120)
        self.bot = DiscordBot.get_instance()
        self.information = {}
        self.author_snowflake = int(interaction.user.id)
        self.interaction = interaction

    async def interaction_check(self, interaction):
        return interaction.user.id == self.author_snowflake

    async def setup(self):
        guild_options = await self._build_guild_options()
        channel_options = await self._build_channel_options()
        duration_options = self._build_duration_options()
        infraction_options = self._build_infraction_options()
        self.guild_select.options = guild_options
        self.channel_select.options = channel_options
        self.duration_select.options = duration_options
        self.infraction_select.options = infraction_options

    async def _build_guild_options(self):
        options = []
        default_kwargs = {"member_snowflake": self.author_snowflake}
        developer = await Developer.select(**default_kwargs)
        sysadmin = None
        try:
            sysadmin = await is_sysadmin_wrapper(source=self.interaction)
        except NotSysadmin:
            pass
        can_all_guilds = bool(developer) or bool(sysadmin)
        if can_all_guilds:
            options.append(discord.SelectOption(label="All", value="all"))
        for g in self.bot.guilds:
            options.append(discord.SelectOption(label=g.name, value=str(g.id)))
        return options

    async def _build_channel_options(self):
        available_channel_snowflakes = set()
        default_kwargs = {"guild_snowflake": int(self.interaction.guild.id), "member_snowflake": self.author_snowflake}
        moderators = await Moderator.select(**default_kwargs)
        if moderators:
            available_channel_snowflakes.update(int(m.channel_snowflake) for m in moderators)
        coordinators = await Coordinator.select(**default_kwargs)
        if coordinators:
            available_channel_snowflakes.update(int(c.channel_snowflake) for c in coordinators)
        administrator = await Administrator.select(**default_kwargs)
        is_admin = bool(administrator)
        guild_owner = None
        try:
            guild_owner = await is_guild_owner_wrapper(source=self.interaction)
        except NotGuildOwner:
            pass
        del default_kwargs["guild_snowflake"]
        developer = await Developer.select(**default_kwargs)
        sysadmin = None
        try:
            sysadmin = await is_sysadmin_wrapper(source=self.interaction)
        except NotSysadmin:
            pass
        can_all_channels = is_admin or guild_owner
        options = []
        if can_all_channels:
            options.append(discord.SelectOption(label="All", value="all"))
        guild = self.interaction.guild
        for ch in guild.channels:
            if isinstance(ch, discord.VoiceChannel) and (ch.id in available_channel_snowflakes or can_all_channels):
                options.append(discord.SelectOption(label=ch.name, value=str(ch.id)))
        return options


    def _build_duration_options(self):
        durations = ["1d", "1w", "1m", "1y"]
        return [
            discord.SelectOption(
                label=duration, value=duration
            )
            for duration in durations
        ]

    def _build_infraction_options(self):
        infractions = [Ban, TextMute, VoiceMute]
        return [
            discord.SelectOption(label=infraction.identifier, value=infraction.identifier)
            for infraction in infractions
            if infraction.identifier != "smute"
        ]

    @discord.ui.select(placeholder="Select guild", options=[])
    async def guild_select(self, interaction, select):
        value = select.values[0]
        self.information["guild_snowflake"] = None if value == "all" else int(value)
        self.guild_select.placeholder = "All" if value == "all" else interaction.guild.name
        self.channel_select.options = await self._build_channel_options()
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(placeholder="Select channel", options=[])
    async def channel_select(self, interaction, select):
        value = select.values[0]
        self.information["channel_snowflake"] = None if value == "all" else int(value)
        self.channel_select.placeholder = "All" if value == "all" else interaction.guild.get_channel(int(value)).name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select duration",
        options=[],
    )
    async def duration_select(self, interaction, select):
        duration_name = select.values[0]
        self.duration_select.placeholder = self.information["duration"] = duration_name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select infraction",
        options=[],
    )
    async def infraction_select(self, interaction, select):
        dir_paths = []
        dir_paths.append(Path("/app/vyrtuous/db/infractions"))
        infraction_name = select.values[0]
        for obj in dir_to_classes(dir_paths=dir_paths, parent=Infraction):
            if obj.identifier == infraction_name and obj.identifier != "smute":
                self.information["infraction"] = obj
                self.information['title'] = obj.__name__
        self.infraction_select.placeholder = infraction_name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        column_names = ["created_at", "infraction_type"]
        conditions = []
        values = []
        if self.information.get("channel_snowflake"):
            conditions.append(f"channel_snowflake=${len(values)+1}")
            values.append(self.information["channel_snowflake"])
        if self.information.get("guild_snowflake"):
            conditions.append(f"guild_snowflake=${len(values)+1}")
            values.append(self.information["guild_snowflake"])
        if "infraction" in self.information:
            conditions.append(f"infraction_type=${len(values)+1}")
            values.append(self.information["infraction"].identifier)
        if "duration" in self.information:
            conditions.append(f"created_at >= ${len(values)+1}")
            values.append(datetime.now(timezone.utc) - DurationObject(self.information["duration"]).to_timedelta())
    
        where_clause = f"WHERE {' AND '.join(conditions)}" if conditions else ""
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch(f"SELECT created_at, infraction_type FROM moderation_logs {where_clause}", *values)
        df = pd.DataFrame(rows, columns=column_names)
        df['created_at'] = pd.to_datetime(df['created_at'], utc=True)
        df = df.sort_values('created_at')
        df['cumulative'] = range(1, len(df) + 1)
        x = mdates.date2num(df['created_at'])
        y = df['cumulative'].values
        coef = np.polyfit(x, y, 1)
        fit = np.poly1d(coef)
        y_pred = fit(x)
        r2 = 1 - ((y - y_pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        sns.set_theme(style='darkgrid', context='talk')
        buffer = BytesIO()
        plt.figure(figsize=(12,6))
        sns.lineplot(data=df, x='created_at', y='cumulative')
        plt.plot(df['created_at'], y_pred)
        plt.title(f'Cumulative {self.information["title"]}s Over Time (R² = {r2:.3f})')
        plt.xlabel('Time')
        plt.ylabel(f'Total {self.information["title"]}s')
        ax = plt.gca()
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=6))
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
        plt.xticks(rotation=30)
        plt.tight_layout()
        plt.savefig(buffer, format='png')
        plt.close()
        buffer.seek(0)

        file = discord.File(fp=buffer, filename='infractions_over_time.png')
        await interaction.response.send_message(file=file)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.message.delete()
        self.stop()
