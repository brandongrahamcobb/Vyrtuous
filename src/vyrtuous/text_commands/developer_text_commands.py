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

from datetime import timedelta
from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database import Database
from vyrtuous.db.developer import NotDeveloper
from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES, at_home
from vyrtuous.listing import list_bugs
from vyrtuous.models.multi_converter import MultiConverter
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.snowflake_context import SnowflakeContext
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.statistics import system_monitoring_service
from vyrtuous.utils.tracking import bug_service
from vyrtuous.utils.users import developer_service, sysadmin_service


class DevTextCommands(commands.Cog):

    PERMISSION_LEVEL = "Developer"

    def __init__(self, *, bot: DiscordBot):
        self.__bot = bot

    async def cog_check(self, ctx: commands.Context):
        context = SnowflakeContext(
            channel_snowflake=ctx.channel.id,
            guild_snowflake=ctx.guild.id,
            member_snowflake=ctx.author.id,
        )
        for verify in (
            sysadmin_service.is_sysadmin_wrapper,
            developer_service.is_developer_wrapper,
        ):
            try:
                if await verify(context=context):
                    return True
            except commands.CheckFailure:
                continue
        raise NotDeveloper

    @commands.command(name="backup", help="DB backup.")
    async def backup_text_command(self, ctx: commands.Context):
        tick = Tick(bot=self.__bot, ctx=ctx)
        db = Database(config=self.__bot.config, directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await tick.end(warning=str(e).capitalize())
        return await tick.end(success=discord.File(db.file_name))

    @commands.command(name="stats", help="Lists stats.")
    async def list_stats_text_command(self, ctx: commands.Context):
        tick = Tick(bot=self.__bot, ctx=ctx)
        embed = discord.Embed(title="Statistics")
        cpu_usage = await system_monitoring_service.calculate_cpu_usage()
        embed.add_field(name="CPU %", value=cpu_usage, inline=True)
        with open("/sys/fs/cgroup/memory.current", "r") as file:
            bits = int(file.read())
            memory_usage = round((bits / 1024) / 1024, 0)
            embed.add_field(name="RAM", value=f"{memory_usage} MB", inline=True)
        with open("/proc/uptime", "r") as file:
            content = file.readline()
            fields = content.split()
            uptime_seconds = float(fields[0])
            time = timedelta(seconds=uptime_seconds)
        embed.add_field(name="Uptime", value=f"{str(time)}", inline=True)
        rx_usage = await system_monitoring_service.calculate_rx_usage()
        embed.add_field(name="RX MB", value=f"{rx_usage} MB", inline=True)
        tx_usage = await system_monitoring_service.calculate_tx_usage()
        embed.add_field(name="TX MB", value=f"{tx_usage} MB", inline=True)
        number_of_servers = len(self.__bot.guilds)
        embed.add_field(name="Servers", value=number_of_servers, inline=True)
        return await tick.end(success=embed)

    @commands.command(
        name="bug", help="Resolve or update the notes on an issue by reference"
    )
    async def update_bug_tracking_text_command(
        self,
        ctx: commands.Context,
        reference: str = commands.parameter(
            description="Specify the developer log reference ID."
        ),
        action: str = commands.parameter(
            description="Specify one of: `resolve` or `append` or `overwrite`.",
        ),
        *,
        notes: str = commands.parameter(
            default=None, description="Optionally specify notes."
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        msg = await bug_service.update_bug(
            action=action, notes=notes, reference=reference
        )
        return await tick.end(success=msg)

    @commands.command(name="bugs", help="List issues.")
    async def list_bugs_text_command(
        self,
        ctx: commands.Context,
        target: Union[str, discord.Guild, None] = commands.parameter(
            converter=MultiConverter,
            default=None,
            description="Specify one of: `all`, server ID or UUID.",
        ),
        *,
        scope: str = commands.parameter(
            default=None,
            description="Optionally specify `resolved` or `unresolved`.",
        ),
    ):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if target == "all":
            obj = None
            reference = None
        elif isinstance(target, discord.Guild):
            obj = target
            reference = None
        else:
            reference = target
            obj = ctx.guild
        is_at_home = at_home(source=ctx)
        pages = await list_bugs.build_pages(
            is_at_home=is_at_home, obj=obj, reference=reference, scope=scope
        )
        return await tick.end(success=pages)

    @commands.command(name="cogs", help="Lists cogs.")
    async def list_cogs_text_command(self, ctx: commands.Context):
        tick = Tick(bot=self.__bot, ctx=ctx)
        if ctx.guild is None:
            return await tick.end(warning="This command must be used in a server.")
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
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        tick = Tick(bot=self.__bot, ctx=ctx)
        try:
            await self.__bot.load_extension(module)
        except commands.ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully loaded {module}.")

    @commands.command(name="ping", help="Ping me!")
    async def ping_text_command(self, ctx: commands.Context):
        tick = Tick(bot=self.__bot, ctx=ctx)
        return await tick.end(success="Pong!")

    @commands.command(
        name="reload", help="Reloads a cog by name `vyrtuous.cog.<cog_name>`."
    )
    async def reload_text_command(self, ctx: commands.Context, *, module: str):
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
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        tick = Tick(bot=self.__bot, ctx=ctx)
        try:
            await self.__bot.unload_extension(module)
        except commands.ExtensionError as e:
            return await tick.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await tick.end(success=f"Successfully unloaded {module}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(DevTextCommands(bot=bot))
