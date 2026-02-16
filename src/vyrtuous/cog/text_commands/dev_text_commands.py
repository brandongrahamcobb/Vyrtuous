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

from typing import Any, Coroutine, Union
from uuid import UUID

import discord
from discord.ext import commands

from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bug.bug_service import BugService
from vyrtuous.database import Database
from vyrtuous.developer.developer import Developer
from vyrtuous.developer.developer_service import DeveloperService, NotDeveloper
from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES
from vyrtuous.sysadmin.sysadmin_service import SysadminService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.discord_object_service import DiscordObjectService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.home import at_home
from vyrtuous.utils.logger import logger
from vyrtuous.utils.message_service import MessageService
from vyrtuous.utils.state_service import StateService


class DevTextCommands(commands.Cog):
    ROLE = Developer

    def __init__(self, bot: DiscordBot):
        self.__bot = bot
        self.message_service = MessageService(self.__bot)
        self.__author_service = AuthorService()
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__bug_service = BugService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__sysadmin_service = SysadminService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
        )
        self.__developer_service = DeveloperService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            emoji=self.__emoji,
        )
        self.__discord_object_service = DiscordObjectService()

    async def cog_check(self, ctx) -> Coroutine[Any, Any, bool]:
        async def predicate(
            source: Union[commands.Context, discord.Interaction, discord.Message],
        ):
            for verify in (
                self.__sysadmin_service.is_sysadmin_wrapper,
                self.__developer_service.is_developer_wrapper,
            ):
                try:
                    if await verify(source):
                        return True
                except commands.CheckFailure:
                    continue
            raise NotDeveloper

        predicate._permission_level = "Developer"
        return await predicate(ctx)

    @commands.command(name="backup", help="DB backup.")
    async def backup_text_command(self, ctx: commands.Context):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        db = Database(config=self.__bot.config, directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await state.end(warning=str(e).capitalize())
        return await state.end(success=discord.File(db.file_name))

    @commands.command(name="cogs", help="Lists cogs.")
    async def list_cogs_text_command(self, ctx: commands.Context):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{self.__emoji.get_random_emoji()} Cogs for {ctx.guild.me.name}",
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
        return await state.end(success=embed)

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
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        msg = self.__bug_service.interact_with_bug(
            action=action, notes=notes, reference=reference
        )
        return await state.end(success=f"{self.__emoji.get_random_emoji()} {msg}.")

    @commands.command(name="bugs", help="List issues.")
    async def list_bugs_text_command(
        self,
        ctx: commands.Context,
        target: Union[str, discord.Guild, None] = commands.parameter(
            converter=commands.GuildConverter,
            default=None,
            description="Specify one of: `all`, server ID or UUID.",
        ),
        *,
        scope: str = commands.parameter(
            default=None,
            description="Optionally specify `resolved` or `unresolved`.",
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        obj = target or ctx.guild
        is_at_home = at_home(source=ctx)
        try:
            target_uuid = UUID(str(target))
            where_kwargs = {"id": target_uuid}
        except Exception as e:
            logger.warning(str(e).capitalize())
            object_dict = await self.__discord_object_service.to_dict(obj=obj)
            where_kwargs = object_dict.get("columns", None)
        pages = await self.__bug_service.build_pages(
            scope=scope, where_kwargs=where_kwargs, is_at_home=is_at_home
        )
        return await state.end(success=pages)

    @commands.command(
        name="load", help="Loads a cog by name 'vyrtuous.cog.<cog_name>.'"
    )
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        try:
            await self.__bot.load_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully loaded {module}.")

    @commands.command(name="ow", help="Overwrite stats.")
    async def list_overwrites_text_command(
        self,
        ctx: commands.Context,
        channel: discord.abc.GuildChannel = commands.parameter(
            converter=commands.VoiceChannelConverter,
            default=None,
            description="Specify the ID or mention.",
        ),
    ):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        obj = channel or ctx.channel
        object_dict = self.__discord_object_service.to_dict(obj=obj)
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
    async def ping_text_command(self, ctx: commands.Context):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        return await state.end(success="Pong!")

    @commands.command(
        name="reload", help="Reloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    async def reload_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        try:
            await self.__bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully reloaded {module}.")

    @commands.command(
        name="unload", help="Unloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(
            author_service=self.__author_service,
            bot=self.__bot,
            bug_service=self.__bug_service,
            ctx=ctx,
            developer_service=self.__developer_service,
            emoji=self.__emoji,
        )
        try:
            await self.__bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully unloaded {module}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(DevTextCommands(bot))
