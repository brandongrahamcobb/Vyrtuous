"""cancel_confirm.py The purpose of this program is to provide an embed utility with cancellation and confirmation buttons.

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

import discord

from vyrtuous.bot.discord_bot import DiscordBot


class VerifyView(discord.ui.View):

    def __init__(
        self,
        category,
        author_snowflake,
        guild_snowflake,
        channel_snowflake=None,
        member_snowflake=None,
        timeout=60,
    ):
        super().__init__(timeout=timeout)
        match category:
            case "alias":
                self.action = "Deletes all aliases."
            case "admin":
                self.action = "Deletes all administrators."
            case "arole":
                self.action = "Deletes all administrator roles."
            case "all":
                self.action = "Deletes all of the above: administrators, administrator roles, aliases, bans, coords, devs, flags, mods, stages, temporary rooms, text-mutes, vegans, voice-mutes and video rooms."
            case "ban":
                self.action = "Deletes all bans."
            case "coord":
                self.action = "Deletes all coordinators."
            case "dev":
                self.action = "Deletes all developers."
            case "flag":
                self.action = "Deletes all flags."
            case "mod":
                self.action = "Deletes all moderators."
            case "stage":
                self.action = "Deletes all stages."
            case "temp":
                self.action = "Deletes all temporary rooms."
            case "tmute":
                self.action = "Deletes all text-mutes."
            case "vegan":
                self.action = "Deletes all vegans."
            case "vmute":
                self.action = "Deletes all voice-mutes."
            case "vr":
                self.action = "Deletes all video rooms."
            case _:
                raise ValueError("Invalid action type specified for confirmation view.")
        self.author_snowflake = author_snowflake
        self.bot = DiscordBot.get_instance()
        self.guild_snowflake = guild_snowflake
        self.channel_snowflake = channel_snowflake
        self.member_snowflake = member_snowflake
        self.result = None
        self.target = self._resolve_target()

    def _resolve_target(self):
        guild = self.bot.get_guild(self.guild_snowflake)
        target = guild.get_member(self.member_snowflake)
        if target:
            return target
        return guild.get_channel(self.channel_snowflake)

    async def interaction_check(self, interaction):
        return interaction.user.id == self.author_snowflake

    @discord.ui.button(label="Confirm", style=discord.ButtonStyle.green)
    async def confirm(self, interaction, button):
        self.result = True
        await interaction.message.delete()
        self.stop()

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        self.result = False
        await interaction.message.delete()
        self.stop()

    def build_embed(self, target):
        mention = target.mention if target else "Unknown"
        embed = discord.Embed(
            title="\U000026a0\U0000fe0f Clear Command Confirmation",
            description=f"**Action:** {self.action}\n**Target:** {mention}",
            color=discord.Color.orange(),
        )
        embed.set_footer(text="Please confirm or cancel this action.")
        return embed
