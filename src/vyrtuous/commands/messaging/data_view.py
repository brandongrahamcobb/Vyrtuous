"""!/bin/python3
data_view.py The purpose of this program is to provide the data view utility class.

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
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.infractions.ban.ban import Ban
from vyrtuous.db.infractions.tmute.text_mute import TextMute
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute
from vyrtuous.utils.dir_to_classes import dir_to_classes


class DataView(discord.ui.View):
    def __init__(self, interaction: discord.Interaction):
        super().__init__(timeout=120)
        self.bot = DiscordBot.get_instance()
        self.information = {}
        self.author_snowflake = int(interaction.user.id)
        self.interaction = interaction

    async def interaction_check(self, interaction):
        return interaction.user.id == self.author_snowflake

    async def setup(self):
        channel_options, guild_options = await self._build_discord_options()
        duration_options = self._build_duration_options()
        infraction_options = self._build_infraction_options()
        self.guild_select.options = guild_options
        self.channel_select.options = channel_options
        self.duration_select.options = duration_options
        self.infraction_select.options = infraction_options

    async def _build_discord_options(self):
        channel_options = []
        guild_options = []
        available_channels, available_guilds = await PermissionService.can_list(
            source=self.interaction
        )
        if "all" in available_channels:
            channel_list = available_channels["all"]
            channel_options.append(discord.SelectOption(label="All", value="all"))
        else:
            channel_list = [
                c
                for ch_list in available_channels.values()
                for c in ch_list
                if c != "all"
            ]
        if "all" in available_guilds:
            guild_list = available_guilds["all"]
            guild_options.append(discord.SelectOption(label="All", value="all"))
        else:
            guild_list = [g for g in available_guilds.values() if g != "all"]
        channel_options.extend(
            discord.SelectOption(label=c.name, value=str(c.id)) for c in channel_list
        )
        guild_options.extend(
            discord.SelectOption(label=g.name, value=str(g.id)) for g in guild_list
        )
        return channel_options, guild_options

    def _build_duration_options(self):
        durations = ["1d", "1w", "1m", "1y", "5y"]
        return [
            discord.SelectOption(label=duration, value=duration)
            for duration in durations
        ]

    def _build_infraction_options(self):
        infractions = [Ban, TextMute, VoiceMute]
        return [
            discord.SelectOption(
                label=infraction.identifier, value=infraction.identifier
            )
            for infraction in infractions
            if infraction.identifier != "smute"
        ]

    @discord.ui.select(placeholder="Select guild", options=[])
    async def guild_select(self, interaction, select):
        value = select.values[0]
        self.information["guild_snowflake"] = None if value == "all" else int(value)
        self.guild_select.placeholder = (
            "All" if value == "all" else interaction.guild.name
        )
        await interaction.response.defer()
        channel_options, _ = await self._build_discord_options()
        if value != "all":
            guild_id = int(value)
            channel_options = [
                discord.SelectOption(label=c.name, value=str(c.id))
                for c in (await interaction.client.fetch_guild(guild_id)).channels
                if isinstance(c, discord.abc.GuildChannel)
            ]
        for item in self.children:
            if isinstance(
                item, discord.ui.Select
            ) and item.placeholder.lower().startswith("select channel"):
                item.options = channel_options
        await interaction.edit_original_response(view=self)

    @discord.ui.select(placeholder="Select channel", options=[])
    async def channel_select(self, interaction, select):
        value = select.values[0]
        self.information["channel_snowflake"] = None if value == "all" else int(value)
        self.channel_select.placeholder = (
            "All" if value == "all" else interaction.guild.get_channel(int(value)).name
        )
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
                self.information["title"] = obj.__name__
        self.infraction_select.placeholder = infraction_name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        column_names = ["created_at", "infraction_type"]
        conditions = []
        values = []
        if self.information.get("channel_snowflake"):
            conditions.append(f"channel_snowflake=${len(values) + 1}")
            values.append(self.information["channel_snowflake"])
        if self.information.get("guild_snowflake"):
            conditions.append(f"guild_snowflake=${len(values) + 1}")
            values.append(self.information["guild_snowflake"])
        if "infraction" in self.information:
            conditions.append(f"infraction_type=${len(values) + 1}")
            values.append(self.information["infraction"].identifier)
        if "duration" in self.information:
            conditions.append(f"created_at >= ${len(values) + 1}")
            values.append(
                datetime.now(timezone.utc)
                - DurationObject(self.information["duration"]).to_timedelta()
            )

        where_clause = f"WHERE {' AND '.join(conditions)}" if conditions else ""
        async with self.bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                f"SELECT created_at, infraction_type FROM moderation_logs {where_clause}",
                *values,
            )
        df = pd.DataFrame(rows, columns=column_names)
        df["created_at"] = pd.to_datetime(df["created_at"], utc=True)
        df = df.sort_values("created_at")
        df["cumulative"] = range(1, len(df) + 1)
        x = mdates.date2num(df["created_at"])
        y = df["cumulative"].values
        coef = np.polyfit(x, y, 1)
        fit = np.poly1d(coef)
        y_pred = fit(x)
        r2 = 1 - ((y - y_pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        sns.set_theme(style="darkgrid", context="talk")
        buffer = BytesIO()
        plt.figure(figsize=(12, 6))
        sns.lineplot(data=df, x="created_at", y="cumulative")
        plt.plot(df["created_at"], y_pred)
        plt.title(f"Cumulative {self.information['title']}s Over Time (R² = {r2:.3f})")
        plt.xlabel("Time")
        plt.ylabel(f"Total {self.information['title']}s")
        ax = plt.gca()
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=6))
        ax.xaxis.set_major_formatter(
            mdates.ConciseDateFormatter(ax.xaxis.get_major_locator())
        )
        plt.xticks(rotation=30)
        plt.tight_layout()
        plt.savefig(buffer, format="png")
        plt.close()
        buffer.seek(0)

        file = discord.File(fp=buffer, filename="infractions_over_time.png")
        await interaction.response.send_message(file=file)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.message.delete()
        self.stop()
