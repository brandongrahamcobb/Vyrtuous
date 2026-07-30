"""!/bin/python3
sysadmin_text_commands.py A discord.py cog containing sysadmin commands for the Vyrtuous bot.

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

from typing import Literal, Optional, Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database import Database
from vyrtuous.models.metadata import metadata
from vyrtuous.models.module import Module, ModuleObject
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import ExtensionError
from vyrtuous.utils.messaging.tick import Tick


class DevelopmentTextCommands(commands.Cog):

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    @commands.command(name="backup", help="Request a backup file.")
    @metadata(permission="command.dev.backup")
    async def backup_text_command(self, ctx: commands.Context) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            requested=["command.dev.backup"],
        )
        db = Database(config=self.__bot.config, directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await tick.end(warning=str(e).capitalize())
        return await tick.end(success=discord.File(db.file_name))

    @commands.command(name="load", help="Loads a cog.")
    @metadata(permission="command.dev.load")
    async def load_text_command(
        self,
        ctx: commands.Context,
        module: ModuleObject = commands.parameter(
            converter=Module, default=None, description="Specify the cog name."
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            requested=["command.dev.load"],
        )
        try:
            await self.__bot.load_extension(module.module)
        except ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully loaded {module}.")

    @commands.command(name="reload", help="Reloads a cog.")
    @metadata(permission="command.dev.reload")
    async def reload_text_command(
        self,
        ctx: commands.Context,
        module: ModuleObject = commands.parameter(
            converter=Module, default=None, description="Specify the cog name."
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            requested=["command.dev.reload"],
        )
        try:
            await self.__bot.reload_extension(module.module)
        except ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully reloaded {module}.")

    @commands.command(name="sync", help="Sync app commands.")
    @metadata(permission="command.dev.sync")
    async def sync_text_command(
        self,
        ctx: commands.Context,
        spec: Optional[Literal["~", "*", "^"]] = None,
        *,
        guilds: Union[commands.Greedy[discord.Object], None] = None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            requested=["command.dev.sync"],
        )
        synced = []
        if not guilds:
            if spec == "~":
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == "*":
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == "^":
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
            else:
                synced = await ctx.bot.tree.sync()
            try:
                if spec is None:
                    msg = f"Synced {len(synced)} commands globally."
                else:
                    msg = f"Synced {len(synced)} commands to the current server."
                return await tick.end(success=msg)
            except Exception as e:
                return await tick.end(warning=str(e).capitalize())
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        return await tick.end(success=f"Synced the tree to {ret}/{len(guilds)}.")

    @commands.command(
        name="unload",
        help="Unloads a cog.",
        extras={"permission": "command.dev.unload"},
    )
    @metadata(permission="command.dev.unload")
    async def unload_text_command(
        self,
        ctx: commands.Context,
        module: ModuleObject = commands.parameter(
            converter=Module, default=None, description="Specify a cog."
        ),
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
        if ctx.channel is None:
            return await tick.end(
                warning="This command must be used in a server channel."
            )
        bot: DiscordBot = DiscordBot.get_instance()
        permission_state: PermissionState = bot.registry.get(PermissionState)
        await permission_service.has_permissions(
            permission_state=permission_state,
            member_snowflake=ctx.author.id,
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            requested=["command.dev.unload"],
        )
        try:
            await self.__bot.unload_extension(module.module)
        except ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully unloaded {module}.")

    # @commands.command(name="upload", help="Create the upload document.")
    # async def uploads_text_command(
    #     self,
    #     ctx: commands.Context,
    # ) -> discord.Message:
    #     tick = Tick(bot=self.__bot, ctx=ctx)
    #     await upload_service.build_latex_document()
    #     return await tick.end(success="Success!")


async def setup(bot: DiscordBot):
    await bot.add_cog(DevelopmentTextCommands(bot=bot))
