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

from io import BytesIO

import discord
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class DataView(discord.ui.View):
    def __init__(
        self,
        *,
        bot=None,
        ban_service=None,
        duration_service=None,
        flag_service=None,
        interaction: discord.Interaction = None,
        moderator_service=None,
        text_mute_service=None,
        voice_mute_service=None,
        state=None,
    ):
        super().__init__(timeout=120)
        self.__bot = bot
        self.__information = {}
        self.__available_channels = {}
        self.__author_snowflake = int(interaction.user.id)
        self.__interaction = interaction
        self.__moderator_service = moderator_service
        self.__ban_service = ban_service
        self.__voice_mute_service = voice_mute_service
        self.__flag_service = flag_service
        self.__text_mute_service = text_mute_service
        self.__infractions = [
            self.__ban_service.MODEL,
            self.__flag_service.MODEL,
            self.__text_mute_service.MODEL,
            self.__voice_mute_service.MODEL,
        ]
        self.__state = state
        self.__duration_service = duration_service

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__author_snowflake

    async def setup(self):
        self.channel_select.options = [
            discord.SelectOption(
                label="Select a guild first",
                value="none",
                default=True,
                description="Please select a guild to see channels",
                emoji=None,
            )
        ]
        duration_options = self._build_duration_options()
        guild_options = await self._build_guild_options()
        infraction_options = self._build_infraction_options()
        self.guild_select.options = guild_options
        self.duration_select.options = duration_options
        self.infraction_select.options = infraction_options

    def _build_duration_options(self):
        durations = ["1d", "1w", "1m", "1y", "5y"]
        return [
            discord.SelectOption(label=duration, value=duration)
            for duration in durations
        ]

    async def _build_guild_options(self):
        guild_options = []
        available_channels, available_guilds = await self.__moderator_service.can_list(
            source=self.__interaction
        )
        self.__available_channels = available_channels
        if "all" in available_guilds:
            guild_list = available_guilds["all"]
            guild_options.append(discord.SelectOption(label="All", value="all"))
        else:
            guild_list = [g for g in available_guilds.values() if g != "all"]
        guild_options.extend(
            discord.SelectOption(label=g.name, value=str(g.id)) for g in guild_list
        )
        return guild_options

    def _build_infraction_options(self):
        return [
            discord.SelectOption(
                label=infraction.identifier, value=infraction.identifier
            )
            for infraction in self.__infractions
            if infraction.identifier != "smute"
        ]

    @discord.ui.select(placeholder="Select guild", options=[])
    async def guild_select(self, interaction, select):
        value = select.values[0]
        self.__information["guild_snowflake"] = None if value == "all" else int(value)
        self.guild_select.placeholder = (
            "All" if value == "all" else interaction.guild.name
        )
        await interaction.response.defer()
        channel_options = []
        if value == "all":
            top_channels = []
            for guild_id, ch_list in self.__available_channels.items():
                if guild_id == "all":
                    continue
                guild = interaction.client.get_guild(int(guild_id))
                if guild:
                    channels = [
                        c for c in guild.channels if isinstance(c, discord.VoiceChannel)
                    ]
                    channels.sort(
                        key=lambda c: getattr(c, "member_count", 0), reverse=True
                    )
                    top_channels.extend(channels[:25])
            top_channels = list({c.id: c for c in top_channels}.values())
            channel_options = [
                discord.SelectOption(label=c.name, value=str(c.id))
                for c in top_channels
            ]
            channel_options.insert(0, discord.SelectOption(label="All", value="all"))
        if value != "all":
            guild_obj = interaction.client.get_guild(int(value))
            channels = [
                c for c in guild_obj.channels if isinstance(c, discord.VoiceChannel)
            ]
            channels.sort(key=lambda c: getattr(c, "member_count", 0), reverse=True)
            top_channels = channels[:25]
            top_ids = {c.id for c in top_channels}
            filtered_channels = []
            seen = set()
            for ch_list in self.__available_channels.values():
                for ch in ch_list:
                    if ch != "all" and ch.id in top_ids and ch.id not in seen:
                        seen.add(ch.id)
                        filtered_channels.append(ch)
            channel_options = [
                discord.SelectOption(label=c.name, value=str(c.id))
                for c in filtered_channels
            ]
            channel_options.insert(0, discord.SelectOption(label="All", value="all"))
        for item in self.children:
            if (
                isinstance(item, discord.ui.Select)
                and item.placeholder
                and item.placeholder.lower().startswith("select channel")
            ):
                item.options = channel_options
        await interaction.edit_original_response(view=self)

    @discord.ui.select(placeholder="Select channel", options=[])
    async def channel_select(self, interaction, select):
        value = select.values[0]
        self.__information["channel_snowflake"] = None if value == "all" else int(value)
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
        self.duration_select.placeholder = self.__information["duration"] = (
            duration_name
        )
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.select(
        placeholder="Select infraction",
        options=[],
    )
    async def infraction_select(self, interaction, select):
        infraction_name = select.values[0]
        for infraction in self.__infractions:
            if infraction.identifier == infraction_name:
                self.__information["infraction"] = infraction
                self.__information["title"] = infraction.__name__
        self.infraction_select.placeholder = infraction_name
        await interaction.response.defer()
        await interaction.edit_original_response(view=self)

    @discord.ui.button(label="Submit", style=discord.ButtonStyle.green)
    async def submit(self, interaction, button):
        column_names = ["created_at", "infraction_type"]
        conditions = []
        values = []
        if self.__information.get("channel_snowflake"):
            conditions.append(f"channel_snowflake=${len(values) + 1}")
            values.append(self.__information["channel_snowflake"])
        if self.__information.get("guild_snowflake"):
            conditions.append(f"guild_snowflake=${len(values) + 1}")
            values.append(self.__information["guild_snowflake"])
        if "infraction" in self.__information:
            conditions.append(f"infraction_type=${len(values) + 1}")
            values.append(self.__information["infraction"].identifier)
        if "duration" in self.__information:
            conditions.append(f"created_at >= ${len(values) + 1}")
            duration = self.__duration_service.parse(self.__information["duration"])
            values.append(self.__duration_service.to_expires_in(duration=duration))
        where_clause = f"WHERE {' AND '.join(conditions)}" if conditions else ""
        async with self.__bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                f"SELECT created_at, infraction_type FROM moderation_logs {where_clause}",
                *values,
            )
        df = pd.DataFrame(rows, columns=column_names)
        df["created_at"] = pd.to_datetime(df["created_at"], utc=True)
        df = df.sort_values("created_at")
        df = df.groupby("created_at", as_index=False).size()
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
        plt.title(
            f"Cumulative {self.__information['title']}s Over Time (R² = {r2:.3f})"
        )
        plt.xlabel("Time")
        plt.ylabel(f"Total {self.__information['title']}s")
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
        return await self.__state.end(success=file)

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        await interaction.message.delete()
        self.stop()
        return await self.__state.end(success="Cancelled action.")
