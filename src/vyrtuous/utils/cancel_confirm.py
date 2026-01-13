"""cancel_confirm.py The purpose of this program is to provide an embed utility with cancellation and confirmation buttons.

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
"""

from vyrtuous.bot.discord_bot import DiscordBot
import discord


class VerifyView(discord.ui.View):

    def __init__(
        self,
        action_type,
        author_snowflake,
        guild_snowflake,
        channel_snowflake=None,
        member_snowflake=None,
        timeout=60,
    ):
        super().__init__(timeout=timeout)
        match action_type:
            case "alias":
                self.action_type = "Deletes all aliases."
            case "all":
                self.action_type = "Deletes all of the above: aliases, bans, coords, flags, mods, temporary rooms, text-mutes, vegans, voice-mutes and video rooms."
            case "ban":
                self.action_type = "Deletes all bans."
            case "coord":
                self.action_type = "Deletes all coordinators."
            case "flag":
                self.action_type = "Deletes all flags."
            case "mod":
                self.action_type = "Deletes all moderators."
            case "temp":
                self.action_type = "Deletes all temporary rooms."
            case "tmute":
                self.action_type = "Deletes all text-mutes."
            case "vegan":
                self.action_type = "Deletes all vegans."
            case "vmute":
                self.action_type = "Deletes all voice-mutes."
            case "vr":
                self.action_type = "Deletes all video rooms."
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

    def build_embed(self, action_type, target):
        mention = target.mention if target else "Unknown"
        embed = discord.Embed(
            title="\U000026a0\U0000fe0f Clear Command Confirmation",
            description=f"**Action**: {action_type}\n**Target**: {mention}",
            color=discord.Color.orange(),
        )
        embed.set_footer(text="Please confirm or cancel this action.")
        return embed
