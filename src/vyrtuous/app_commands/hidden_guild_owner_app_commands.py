"""!/bin/python3
guild_owner_app_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

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
from discord import app_commands
from discord.ext import commands

from vyrtuous.app_commands.help_app_command import skip_app_command_help_discovery
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import list_heroes
from vyrtuous.models import multi_converter
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import hero_service, moderator_service


class HiddenGuildOwnerAppCommands(commands.Cog):
    PERMISSION_LEVEL = "Guild Owner"

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def interaction_check(self, interaction: discord.Interaction):
        if interaction.guild is None:
            raise commands.CheckFailure(
                "This command must be executed inside a server."
            )
        if interaction.channel is None:
            raise commands.CheckFailure(
                "This command must be executed in a valid channel."
            )
        await moderator_service.check_minimum_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @app_commands.command(name="hero", description="Grant/revoke invincibility.")
    @app_commands.describe(
        member="Specify a member ID/mention.", server="Specify a server ID"
    )
    @skip_app_command_help_discovery()
    async def toggle_invincibility_app_command(
        self, interaction: discord.Interaction, member: str, server: str | None
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        if server == "all":
            guilds = self.__bot.guilds
        elif server is None:
            guilds = [interaction.guild]
        else:
            guild = multi_converter.transform(interaction=interaction, argument=server)
            if not isinstance(guild, discord.Guild):
                return await tick.end(
                    warning="This command must be executed for a valid server."
                )
            guilds = [guild]
        msg = "No hero granted."
        for guild in guilds:
            resolved = guild.get_member(member_snowflake)
            if member_set := self.__bot.registry.get(MemberState).invincible.get(
                guild.id
            ):
                pass
            else:
                member_set = self.__bot.registry.get(MemberState).invincible[
                    guild.id
                ] = set()
            if member_snowflake not in member_set:
                await hero_service.add_invincible_member(guild.id, member_snowflake)
                await hero_service.unrestrict(
                    guild_snowflake=guild.id, member_snowflake=member_snowflake
                )
                resolved = guild.get_member(member_snowflake)
                msg = (
                    f"All moderation events have been forgiven "
                    f"and invincibility has been enabled for {resolved.mention if resolved else member_snowflake}."
                )
            else:
                await hero_service.remove_invincible_member(guild.id, member_snowflake)
                msg = f"Invincibility has been disabled for {resolved.mention if resolved else member_snowflake}."
        return await tick.end(success=msg)

    @app_commands.command(name="heroes", description="List heroes.")
    @skip_app_command_help_discovery()
    async def list_heroes_app_command(
        self,
        interaction: discord.Interaction,
        *,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.guild
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_heroes.build_pages(
            is_at_home=is_at_home,
            obj=obj,
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenGuildOwnerAppCommands(bot=bot))
