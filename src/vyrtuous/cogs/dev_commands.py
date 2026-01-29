"""dev_commands.py A discord.py cog containing developer commands for the Vyrtuous bot.

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

from typing import Literal, Optional
from uuid import UUID

from discord import app_commands
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database import Database
from vyrtuous.db.mgmt.bug import Bug
from vyrtuous.db.roles.developer import developer_predicator
from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES
from vyrtuous.utils.home import at_home
from vyrtuous.utils.logger import logger
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import (
    DiscordObject,
)
from vyrtuous.utils.emojis import get_random_emoji


class DevCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    # DONE
    @app_commands.command(name="backup", description="DB backup.")
    @developer_predicator()
    async def backup_app_command(self, interaction: discord.Interaction):
        await interaction.response.defer(ephemeral=True)
        state = StateService(source=interaction)
        db = Database(config=self.bot.config, directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await state.end(warning=str(e).capitalize())
        return await state.end(success=discord.File(db.file_name))

    # DONE
    @commands.command(name="backup", help="DB backup.")
    @developer_predicator()
    async def backup_text_command(self, ctx: commands.Context):
        state = StateService(source=ctx)
        db = Database(config=self.bot.config, directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await state.end(warning=str(e).capitalize())
        return await state.end(success=discord.File(db.file_name))

    @app_commands.command(name="cogs", description="Lists cogs.")
    @developer_predicator()
    async def list_cogs_app_command(self, interaction: discord.Interaction):
        state = StateService(source=interaction)
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"Cogs for {interaction.guild.me.name}",
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
            embed.add_field(name="No cogs available.", inline=False)
        return await state.end(success=embed)

    @commands.command(name="cogs", help="Lists cogs.")
    @developer_predicator()
    async def list_cogs_text_command(self, ctx: commands.Context):
        state = StateService(source=ctx)
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
            embed.add_field(name="No cogs available.", inline=False)
        return await state.end(success=embed)

    @app_commands.command(
        name="bug", description="Resolve or update the notes on an issue by reference."
    )
    @app_commands.describe(
        reference="Specify the issue reference ID.",
        action="'resolve' or 'append' or 'overwrite'.",
        notes="Optionally specify notes to append or overwrite.",
    )
    @developer_predicator()
    async def update_bug_tracking_app_command(
        self,
        ctx: commands.Context,
        reference: str,
        action: str,
        notes: str = None,
    ):
        state = StateService(source=ctx)

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
            await bug.append(notes)
            detail = "appended to the previous notes."
        elif action and action.lower() == "overwrite":
            where_kwargs = {"id": reference}
            set_kwargs = {"notes": notes}
            await Bug.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
            detail = "overwrote the previous notes."
        return await state.end(
            success=f"\U000026a0\U0000fe0f " f"You successfully {detail}."
        )

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
        notes: Optional[str] = commands.parameter(
            default=None, description="Optionally specify notes."
        ),
    ):
        state = StateService(source=ctx)

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
            await bug.append(notes)
            detail = "appended to the previous notes."
        elif action and action.lower() == "overwrite":
            where_kwargs = {"id": reference}
            set_kwargs = {"notes": notes}
            await Bug.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
            detail = "overwrote the previous notes"
        return await state.end(
            success=f"\U000026a0\U0000fe0f " f"You successfully {detail}."
        )

    @app_commands.command(name="bugs", description="List issues.")
    @app_commands.describe(
        target="Specify one of: `all`, server ID or UUID.",
        filter="Optionally specify `resolved` or `unresolved`.",
    )
    @developer_predicator()
    async def list_bugs_app_command(
        self,
        interaction: discord.Interaction,
        target: str = None,
        filter: str = None,
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        is_at_home = at_home(source=interaction)
        try:
            target_uuid = UUID(str(target))
            kwargs = {"id": target_uuid}
        except Exception as e:
            logger.warning(str(e).capitalize())
            object_dict = await do.determine_from_target(target=target)
            kwargs = object_dict.get("columns", None)
        pages = await Bug.build_pages(
            filter=filter, kwargs=kwargs, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=Bug.PLURAL, pages=pages, state=state)

    @commands.command(name="bugs", help="List issues.")
    @developer_predicator()
    async def list_bugs_text_command(
        self,
        ctx: commands.Context,
        target: str = commands.parameter(
            description="Specify one of: `all`, server ID or UUID.",
        ),
        *,
        filter: Optional[str] = commands.parameter(
            default=None,
            description="Optionally specify `resolved` or `unresolved`.",
        ),
    ):
        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        is_at_home = at_home(source=ctx)
        try:
            target_uuid = UUID(str(target))
            kwargs = {"id": target_uuid}
        except Exception as e:
            logger.warning(str(e).capitalize())
            object_dict = await do.determine_from_target(target=target)
            kwargs = object_dict.get("columns", None)
        pages = await Bug.build_pages(
            filter=filter, kwargs=kwargs, is_at_home=is_at_home
        )
        await StateService.send_pages(plural=Bug.PLURAL, pages=pages, state=state)

    # DONE
    @app_commands.command(
        name="load", description="Loads a cog by name 'vyrtuous.cog.<cog_name>.'"
    )
    @developer_predicator()
    async def load_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = StateService(source=interaction)
        try:
            await interaction.client.load_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully loaded {module}.")

    # DONE
    @commands.command(
        name="load", help="Loads a cog by name 'vyrtuous.cog.<cog_name>.'"
    )
    @developer_predicator()
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(source=ctx)
        try:
            await self.bot.load_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully loaded {module}.")
    
    @app_commands.command(name='ow', description="List overwrites for a channel.")
    @app_commands.describe(channel="Specify the ID or mention.")
    @developer_predicator()
    async def list_overwrites_app_command(
        self,
        interaction: discord.Interaction,
        channel: str
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
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
        embed.add_field(name="Channel", value=object_dict.get("name", None), inline=False)
        embed.add_field(name="Role overwrites", value=str(role_count), inline=False)
        embed.add_field(name="Member overwrites", value=str(member_count), inline=False)
        embed.add_field(name="Total overwrites", value=str(total_count), inline=False)
        return await state.end(success=embed)

    @commands.command(name='ow')
    @developer_predicator()
    async def overwrites(
        self,
        ctx: commands.Context,
        channel: str = commands.parameter(
            default=None,
            description="Specify the ID or mention.",
        )
    ):
        state = StateService(source=ctx)
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
        embed.add_field(name="Channel", value=object_dict.get("name", None), inline=False)
        embed.add_field(name="Role overwrites", value=str(role_count), inline=False)
        embed.add_field(name="Member overwrites", value=str(member_count), inline=False)
        embed.add_field(name="Total overwrites", value=str(total_count), inline=False)
        return await state.end(success=embed)


    # DONE
    @app_commands.command(name="ping", description="Ping me!")
    @developer_predicator()
    async def ping_app_command(self, interaction: discord.Interaction):
        state = StateService(source=interaction)
        return await state.end(success="Pong!")

    # DONE
    @commands.command(name="ping", help="Ping me!")
    @developer_predicator()
    async def ping_text_command(self, ctx: commands.Context):
        state = StateService(source=ctx)
        return await state.end(success="Pong!")

    # DONE
    @app_commands.command(
        name="reload", description="Reloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @app_commands.check(at_home)
    @developer_predicator()
    async def reload_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = StateService(source=interaction)
        try:
            await interaction.client.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully reloaded {module}.")

    # DONE
    @commands.command(
        name="reload", help="Reloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @developer_predicator()
    async def reload_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(source=ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully reloaded {module}.")

    # DONE
    @app_commands.command(name="sync", description="Sync app commands.")
    @developer_predicator()
    async def sync_app_command(
        self,
        interaction: discord.Interaction,
        spec: Optional[Literal["~", "*", "^"]] = None,
    ):
        await interaction.response.defer(ephemeral=True)
        state = StateService(source=interaction)
        guilds = interaction.client.guilds
        synced = []
        if not guilds:
            if spec == "~":
                synced = await interaction.client.tree.sync(guild=interaction.guild)
            elif spec == "*":
                interaction.client.tree.copy_global_to(guild=interaction.guild)
                synced = await interaction.client.tree.sync(guild=interaction.guild)
            elif spec == "^":
                interaction.client.tree.clear_commands(guild=interaction.guild)
                await interaction.client.tree.sync(guild=interaction.guild)
            else:
                synced = await interaction.client.tree.sync()
            try:
                if spec is None:
                    msg = f"Synced {len(synced)} " f"commands globally."
                else:
                    msg = f"Synced {len(synced)} " f"commands to the current server."
                return await state.end(success=msg)
            except Exception as e:
                return await state.end(warning=str(e).capitalize())
        ret = 0
        for guild in guilds:
            try:
                await interaction.client.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        return await state.end(success=f"Synced the tree to {ret}/{len(guilds)}.")

    # DONE
    @commands.command(name="sync", help="Sync app commands.")
    @developer_predicator()
    async def sync_text_command(
        self,
        ctx: commands.Context,
        spec: Optional[Literal["~", "*", "^"]] = None,
        *,
        guilds: commands.Greedy[discord.Object] = None,
    ):
        state = StateService(source=ctx)
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
                    msg = f"Synced {len(synced)} commands to the " f"current server."
                return await state.end(success=msg)
            except Exception as e:
                return await state.end(warning=str(e).capitalize())
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        return await state.end(success=f"Synced the tree to {ret}/{len(guilds)}.")

    # DONE
    @app_commands.command(
        name="unload", description="Unloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @developer_predicator()
    async def unload_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = StateService(source=interaction)
        try:
            await interaction.client.unload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully unloaded {module}.")

    # DONE
    @commands.command(
        name="unload", help="Unloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @developer_predicator()
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(source=ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            return await state.end(
                warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
            )
        return await state.end(success=f"Successfully unloaded {module}.")


async def setup(bot: DiscordBot):
    await bot.add_cog(DevCommands(bot))
