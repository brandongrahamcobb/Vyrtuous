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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.listing import list_developers
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import administrator_role_service, moderator_service


class GuildOwnerAppCommands(commands.Cog):
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

    @app_commands.command(name="admin", description="Toggle administrator role.")
    @app_commands.describe(
        role="Specify a role ID/mention.", guild="Specify a server ID."
    )
    async def toggle_administrator_by_role_app_command(
        self,
        interaction: discord.Interaction,
        role: app_commands.Transform[TargetObject, AppTarget],
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
        if isinstance(role.target, discord.Role):
            role_snowflake = role.target.id
        else:
            return await tick.end(warning="This command must target a valid role.")
        message = await interaction.original_response()
        message_snowflake = message.id
        pages = await administrator_role_service.toggle_administrator_role(
            author_snowflake=interaction.user.id,
            guild_snowflake=guild_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message.channel.id,
            role_snowflake=role_snowflake,
        )
        return await tick.end(success=pages)

    @app_commands.command(name="devs", description="List devs.")
    @app_commands.describe(member="Specify a member ID/mention.")
    async def list_developers_app_command(
        self,
        interaction: discord.Interaction,
        member: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if member is None:
            obj = None
        else:
            obj = member.target
        pages = await list_developers.build_pages(obj=obj)
        return await tick.end(success=pages)

    # @app_commands.command(name="sync", description="Sync app commands.")
    # async def sync_app_command(
    #     self,
    #     interaction: discord.Interaction,
    #     spec: Optional[Literal["~", "*", "^"]] = None,
    #     *,
    #     guilds: Union[commands.Greedy[discord.Object], None] = None,
    # ) -> discord.Message:
    #     tick = Tick(bot=self.__bot, interaction=interaction)
    #     synced = []
    #     if not guilds:
    #         if spec == "~":
    #             synced = await self.__bot.tree.sync(guild=interaction.guild)
    #         elif spec == "*":
    #             if interaction.guild is None:
    #                 return await tick.end(
    #                     warning="This command must be executed in a server."
    #                 )
    #             self.__bot.tree.copy_global_to(guild=interaction.guild)
    #             synced = await self.__bot.tree.sync(guild=interaction.guild)
    #         elif spec == "^":
    #             self.__bot.tree.clear_commands(guild=interaction.guild)
    #             await self.__bot.tree.sync(guild=interaction.guild)
    #         else:
    #             synced = await self.__bot.tree.sync()
    #         try:
    #             if spec is None:
    #                 msg = f"Synced {len(synced)} commands globally."
    #             else:
    #                 msg = f"Synced {len(synced)} commands to the current server."
    #             return await tick.end(success=msg)
    #         except Exception as e:
    #             return await tick.end(warning=str(e).capitalize())
    #     ret = 0
    #     for guild in guilds:
    #         try:
    #             await self.__bot.tree.sync(guild=guild)
    #         except discord.HTTPException:
    #             pass
    #         else:
    #             ret += 1
    #     return await tick.end(success=f"Synced the tree to {ret}/{len(guilds)}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerAppCommands(bot=bot))
