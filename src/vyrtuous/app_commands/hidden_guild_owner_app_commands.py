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
from vyrtuous.models.target import AppTarget, TargetObject
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
        member="Specify a member ID/mention.",
        enabled="Speficy True or False.",
        target="Specify `all` or a server ID.",
    )
    @skip_app_command_help_discovery()
    async def toggle_invincibility_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject, AppTarget],
        enabled: bool,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        bot: DiscordBot = DiscordBot.get_instance()
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guilds = [interaction.guild]
        else:
            if isinstance(target.target, discord.Guild):
                guilds = [target.target]
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if isinstance(member, int):
            member_snowflake = member
            member_name = "Unknown"
            simplified_member = self.__bot.registry.get(MemberState).active.get(
                member_snowflake
            )
            if simplified_member:
                member_name = simplified_member[0]
        elif isinstance(member, discord.Member):
            member_snowflake = member.id
            member_name = member.mention
        else:
            return await tick.end(warning=f"This command must target a valid member.")
        if enabled:
            header = f"All moderation events have been forgiven and invincibility has been enabled for {member_name} in "
        else:
            header = f"Invincibility has been disabled for {member_name} in "
        guild_names = []
        for guild in guilds:
            if enabled:
                await hero_service.add_invincible_member(guild.id, member_snowflake)
                await hero_service.unrestrict(
                    guild_snowflake=guild.id, member_snowflake=member_snowflake
                )
                guild_names.append(guild.name)
            else:
                await hero_service.remove_invincible_member(guild.id, member_snowflake)
                guild_names.append(guild.name)
        msg = header + ", ".join(guild_names) + "."
        return await tick.end(success=msg)

    @app_commands.command(name="heroes", description="List heroes.")
    @app_commands.describe(
        target="Specify a member ID/mention or server ID.", guild="Specify a server ID."
    )
    @skip_app_command_help_discovery()
    async def list_heroes_app_command(
        self,
        interaction: discord.Interaction,
        target: app_commands.Transform[TargetObject | None, AppTarget] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if guild is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            guild_snowflake = interaction.guild.id
        else:
            if isinstance(guild.target, discord.Guild):
                guild_snowflake = guild.target.id
            else:
                return await tick.end(
                    warning="This command must target a valid server."
                )
        if target is None:
            if interaction.guild is None:
                return await tick.end(
                    warning="This command must target a valid server."
                )
            obj = interaction.guild
        else:
            obj = target.target
        is_at_home = at_home(source=interaction)
        pages = await list_heroes.build_pages(
            guild_snowflake=guild_snowflake,
            is_at_home=is_at_home,
            obj=obj,
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenGuildOwnerAppCommands(bot=bot))
