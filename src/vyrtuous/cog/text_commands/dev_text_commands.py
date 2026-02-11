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

from uuid import UUID

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug import Bug
from vyrtuous.bug.bug_service import BugService
from vyrtuous.database import Database
from vyrtuous.developer.developer import Developer
from vyrtuous.developer.developer_service import developer_predicator
from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES
from vyrtuous.utils.discord_object_service import DiscordObject
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.home import at_home
from vyrtuous.utils.logger import logger
from vyrtuous.utils.message_service import MessageService
from vyrtuous.utils.state_service import StateService


class DevTextCommands(commands.Cog):

    ROLE = Developer

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @commands.command(name="backup", help="DB backup.")
    @developer_predicator()
    async def backup_text_command(self, ctx: commands.Context):
        state = StateService(ctx=ctx)
        db = Database(config=self.bot.config, directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await state.end(warning=str(e).capitalize())
        return await state.end(success=discord.File(db.file_name))

    @commands.command(name="cogs", help="Lists cogs.")
    @developer_predicator()
    async def list_cogs_text_command(self, ctx: commands.Context):
        state = StateService(ctx=ctx)
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"Cogs for {ctx.guild.me.name}",
            color=discord.Color.blurple(),
        )
        for representation, cog in zip(
            sorted(DISCORD_COGS), sorted(DISCORD_COGS_CLASSES)
        ):
            if cog in self.bot.cogs:
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
        return await state.end(success=embed)

    @commands.command(
        name="bug", help="Resolve or update the notes on an issue by reference"
    )
    @developer_predicator()
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
        state = StateService(ctx=ctx)

        bug = await Bug.select(id=reference, resolved=False, singular=True)
        if not bug:
            return await state.end(
                warning=f"Unresolved issue not found for reference ({reference})."
            )
        if action and action.lower() == "resolve":
            where_kwargs = {"id": reference}
            set_kwargs = {"resolved": True}
            await Bug.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
            detail = (
                "resolved the issue. The record "
                "will remain in the database for the next 30 days."
            )
        elif action and action.lower() == "append":
            where_kwargs = {"id": reference}
            set_kwargs = {"notes": bug.notes + notes if bug.notes else notes}
            await Bug.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
            detail = "appended to the previous notes."
        elif action and action.lower() == "overwrite":
            where_kwargs = {"id": reference}
            set_kwargs = {"notes": notes}
            await Bug.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
            detail = "overwrote the previous notes"
        return await state.end(
            success=f"\U000026a0\U0000fe0f " f"You successfully {detail}."
        )

    @commands.command(name="bugs", help="List issues.")
    @developer_predicator()
    async def list_bugs_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: `all`, server ID or UUID.",
        ),
        *,
        filter: str = commands.parameter(
            default=None,
            description="Optionally specify `resolved` or `unresolved`.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        try:
            target_uuid = UUID(str(target))
            where_kwargs = {"id": target_uuid}
        except Exception as e:
            logger.warning(str(e).capitalize())
            object_dict = await do.determine_from_target(target=target)
            where_kwargs = object_dict.get("columns", None)
        pages = await BugService.build_pages(
            filter=filter, where_kwargs=where_kwargs, is_at_home=is_at_home
        )
        await StateService.send_pages(title="Bugs", pages=pages, state=state)

    @commands.command(
        name="load", help="Loads a cog by name 'vyrtuous.cog.<cog_name>.'"
    )
    @developer_predicator()
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(ctx=ctx)
        try:
            await self.bot.load_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully loaded {module}.")

    @commands.command(name="ow", help="Overwrite stats.")
    @developer_predicator()
    async def list_overwrites_text_command(
        self,
        ctx: commands.Context,
        channel: str = commands.parameter(
            default=None,
            description="Specify the ID or mention.",
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        object_dict = await do.determine_from_target(target=channel)
        member_count, role_count, total_count = 0, 0, 0
        for target, overwrite in object_dict.get("object", None).overwrites.items():
            if any(v is not None for v in overwrite):
                total_count += 1
                if isinstance(target, discord.Member):
                    member_count += 1
                else:
                    role_count += 1
        embed = discord.Embed(title="Channel Overwrites", color=discord.Color.blue())
        embed.add_field(
            name="Channel", value=object_dict.get("name", None), inline=False
        )
        embed.add_field(name="Role overwrites", value=str(role_count), inline=False)
        embed.add_field(name="Member overwrites", value=str(member_count), inline=False)
        embed.add_field(name="Total overwrites", value=str(total_count), inline=False)
        return await state.end(success=embed)

    @commands.command(name="ping", help="Ping me!")
    @developer_predicator()
    async def ping_text_command(self, ctx: commands.Context):
        state = StateService(ctx=ctx)
        return await state.end(success="Pong!")

    @commands.command(
        name="reload", help="Reloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @developer_predicator()
    async def reload_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(ctx=ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully reloaded {module}.")

    @commands.command(
        name="unload", help="Unloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @developer_predicator()
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(ctx=ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully unloaded {module}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(DevTextCommands(bot))
