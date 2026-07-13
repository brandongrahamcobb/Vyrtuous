"""!/bin/python3
dev_text_commands.py A discord.py cog containing developer commands for the Vyrtuous bot.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database import Database
from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES
from vyrtuous.text_commands.help_text_command import skip_text_command_help_discovery
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.users import moderator_service


class HiddenDeveloperTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Developer"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx: commands.Context) -> bool:
        if ctx.guild is None:
            raise commands.CheckFailure("This command must be executed inside a server.")
        await moderator_service.check_minimum_role(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @commands.command(name="backup", help="DB backup.")
    @skip_text_command_help_discovery()
    async def backup_text_command(self, ctx: commands.Context) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        db = Database(config=self.__bot.config, directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await tick.end(warning=str(e).capitalize())
        return await tick.end(success=discord.File(db.file_name))

    @commands.command(name="cogs", help="Lists cogs.")
    @skip_text_command_help_discovery()
    async def list_cogs_text_command(self, ctx: commands.Context) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} Cogs for {ctx.guild.me.name}",
            color=discord.Color.blurple(),
        )
        for representation, cog in zip(
            sorted(DISCORD_COGS), sorted(DISCORD_COGS_CLASSES)
        ):
            if cog in self.__bot.cogs:
                loaded.append(representation)
            else:
                not_loaded.append(representation)
        if loaded:
            embed.add_field(name="Loaded", value="\n".join(loaded), inline=False)
        if not_loaded:
            embed.add_field(
                name="Not Loaded", value="\n".join(not_loaded), inline=False
            )
        if not loaded and not not_loaded:
            embed.add_field(name="No cogs available.", value=None, inline=False)
        return await tick.end(success=embed)

    @commands.command(
        name="load", help="Loads a cog by name `vyrtuous.cog.<cog_name>.`"
    )
    @skip_text_command_help_discovery()
    async def load_text_command(
        self, ctx: commands.Context, *, module: str
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        try:
            await self.__bot.load_extension(module)
        except commands.ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully loaded {module}.")

    @commands.command(
        name="reload", help="Reloads a cog by name `vyrtuous.cog.<cog_name>`."
    )
    @skip_text_command_help_discovery()
    async def reload_text_command(
        self, ctx: commands.Context, *, module: str
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        try:
            await self.__bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully reloaded {module}.")

    @commands.command(
        name="unload", help="Unloads a cog by name `vyrtuous.cog.<cog_name>`."
    )
    @skip_text_command_help_discovery()
    async def unload_text_command(
        self, ctx: commands.Context, *, module: str
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, ctx=ctx)
        try:
            await self.__bot.unload_extension(module)
        except commands.ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully unloaded {module}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(HiddenDeveloperTextCommands(bot=bot))
