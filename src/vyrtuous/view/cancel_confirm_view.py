"""!/bin/python3
cancel_confirm.py The purpose of this program is to provide an embed utility with cancellation and confirmation buttons.

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
from discord.ext import commands


class VerifyView(discord.ui.View):
    def __init__(
        self,
        *,
        category=None,
        author_snowflake=None,
        obj=None,
        timeout=60,
    ):
        super().__init__(timeout=timeout)
        match category:
            case "alias":
                self.__action = "Deletes all aliases."
            case "all":
                self.__action = "Deletes all of the above: administrators, administrator roles, aliases, automute channels, bans, coords, devs, flags, mods, text-mutes, vegans, voice-mutes and video channels."
            case "arole":
                self.__action = "Deletes all administrator roles."
            case "automute":
                self.__action = "Deletes all automute channels."
            case "ban":
                self.__action = "Deletes all bans."
            case "coord":
                self.__action = "Deletes all coordinators."
            case "dev":
                self.__action = "Deletes all developers."
            case "flag":
                self.__action = "Deletes all flags."
            case "mod":
                self.__action = "Deletes all moderators."
            case "tmute":
                self.__action = "Deletes all text-mutes."
            case "vegan":
                self.__action = "Deletes all vegans."
            case "vmute":
                self.__action = "Deletes all voice-mutes."
            case "voice_channel":
                self.__action = "Deletes all video channels."
            case _:
                raise commands.CheckFailure(
                    f"Invalid category specified for clear command: {category}."
                )
        self.__author_snowflake = author_snowflake
        if isinstance(obj, discord.Guild):
            self._guild_snowflake = obj.id
        elif isinstance(obj, discord.abc.GuildChannel):
            self._channel_snowflake = obj.id
        elif isinstance(obj, discord.Member):
            self._member_snowflake = obj.id
        self._result = None
        if isinstance(obj, discord.Guild):
            self.__mention = obj.name
        elif obj is not None and not isinstance(obj, str):
            self.__mention = obj.mention
        else:
            self.__mention = "All"

    async def interaction_check(self, interaction):
        return interaction.user.id == self.__author_snowflake

    @discord.ui.button(label="Confirm", style=discord.ButtonStyle.green)
    async def confirm(self, interaction, button):
        self._result = True
        await interaction.message.delete()
        self.stop()

    @discord.ui.button(label="Cancel", style=discord.ButtonStyle.red)
    async def cancel(self, interaction, button):
        self._result = False
        await interaction.message.delete()
        self.stop()

    def build_embed(self):
        embed = discord.Embed(
            title="\U000026a0\U0000fe0f Clear Command Confirmation",
            description=f"**Action:** {self.__action}\n**Target:** {self.__mention}",
            color=discord.Color.orange(),
        )
        embed.set_footer(text="Please confirm or cancel this action.")
        return embed

    @property
    def channel_snowflake(self):
        return self._channel_snowflake

    @property
    def guild_snowflake(self):
        return self._guild_snowflake

    @property
    def member_snowflake(self):
        return self._member_snowflake

    @property
    def result(self):
        return self._result
