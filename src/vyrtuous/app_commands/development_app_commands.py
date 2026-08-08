"""!/bin/python3
hidden_dev_app_commands.py A discord.py cog containing developer commands for the Vyrtuous bot.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from typing import Literal, Optional

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database import Database
from vyrtuous.models.metadata import metadata
from vyrtuous.models.module import AppModule, ModuleObject
from vyrtuous.models.target import AppTarget, TargetObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure, ExtensionError
from vyrtuous.utils.messaging.tick import Tick


class DevelopmentAppCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def cog_app_command_error(self, interaction, error):
        self.__bot.logger.info(str(error))
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.end(error=str(error))

    @metadata(permission="command.dev.backup")
    @app_commands.command(name="backup", description="Request a backup file.")
    async def backup_app_command(
        self, interaction: discord.Interaction
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel.",
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            requested=["command.dev.backup"],
        )
        db = Database(config=self.__bot.config, directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await tick.end(warning=str(e).capitalize())
        return await tick.end(success=discord.File(db.file_name))

    @metadata(permission="command.dev.load")
    @app_commands.command(name="load", description="Loads a cog.")
    @app_commands.describe(module="Specify the cog name.")
    async def load_app_command(
        self,
        interaction: discord.Interaction,
        module: app_commands.Transform[ModuleObject, AppModule],
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            requested=["command.dev.load"],
        )
        try:
            await self.__bot.load_extension(module.module)
        except ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully loaded {module.module}.")

    @metadata(permission="command.dev.reload")
    @app_commands.command(name="reload", description="Reloads a cog.")
    @app_commands.describe(module="Specify the cog name.")
    async def reload_app_command(
        self,
        interaction: discord.Interaction,
        module: app_commands.Transform[ModuleObject, AppModule],
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            requested=["command.dev.reload"],
        )
        try:
            await self.__bot.reload_extension(module.module)
        except ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully reloaded {module.module}.")

    @metadata(permission="command.dev.sync")
    @app_commands.command(name="sync", description="Sync app commands.")
    @app_commands.describe(
        spec="Specify directly to the guild (~), global to guild (*), clear and sync (^) and global sync (None).",
        guild="Specify which guild to sync.",
    )
    async def sync_app_command(
        self,
        interaction: discord.Interaction,
        spec: Optional[Literal["~", "*", "^"]] = None,
        guild: app_commands.Transform[TargetObject | None, AppTarget] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            requested=["command.dev.sync"],
        )
        ret = 0
        synced = []
        if not guild:
            if spec == "~":
                synced = await self.__bot.tree.sync(guild=interaction.guild)
            elif spec == "*":
                if interaction.guild is None:
                    return await tick.end(
                        warning="This command must be executed in a server."
                    )
                self.__bot.tree.copy_global_to(guild=interaction.guild)
                synced = await self.__bot.tree.sync(guild=interaction.guild)
            elif spec == "^":
                self.__bot.tree.clear_commands(guild=interaction.guild)
                await self.__bot.tree.sync(guild=interaction.guild)
            else:
                synced = await self.__bot.tree.sync()
            try:
                if spec is None:
                    msg = f"Synced {len(synced)} commands globally."
                else:
                    msg = f"Synced {len(synced)} commands to the current server."
                return await tick.end(success=msg)
            except Exception as e:
                return await tick.end(warning=str(e).capitalize())
        else:
            if isinstance(guild.target, discord.Guild):
                guild_obj = guild.target
            else:
                raise CheckFailure("This command must target a valid server.")
            try:
                await self.__bot.tree.sync(guild=guild_obj)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        return await tick.end(success=f"Synced the tree to {ret}.")

    # @app_commands.command(name="upload", description="Create the upload document.")
    # async def uploads_app_command(self, interaction: discord.Interaction) -> None:
    #     return await upload_service.build_latex_document()

    @metadata(permission="command.dev.unload")
    @app_commands.command(name="unload", description="Unloads a cog.")
    @app_commands.describe(module="Specify the cog name.")
    async def unload_app_command(
        self,
        interaction: discord.Interaction,
        module: app_commands.Transform[ModuleObject, AppModule],
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        await tick.defer()
        if interaction.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if interaction.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=interaction.user.id,
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            requested=["command.dev.unload"],
        )
        try:
            await self.__bot.unload_extension(module.module)
        except ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully unloaded {module.module}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(DevelopmentAppCommands(bot=bot))
