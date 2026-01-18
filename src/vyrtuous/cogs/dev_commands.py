"""dev_commands.py A discord.py cog containing developer commands for the Vyrtuous bot.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
from vyrtuous.database.database import Database
from vyrtuous.database.logs.developer_log import DeveloperLog
from vyrtuous.inc.helpers import DISCORD_COGS_CLASSES
from vyrtuous.service.check_service import at_home, developer_predicator
from vyrtuous.service.logging_service import logger
from vyrtuous.service.messaging.message_service import MessageService
from vyrtuous.service.messaging.state_service import StateService
from vyrtuous.service.scope_service import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    resolve_objects,
    generate_skipped_guilds,
    generate_skipped_messages,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji


class DevCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot, self.bot.db_pool)

    # DONE
    @app_commands.command(name="backup", description="DB backup.")
    @developer_predicator()
    async def app_backup(self, interaction: discord.Interaction):
        await interaction.response.defer(ephemeral=True)
        state = StateService(interaction)
        db = Database(directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except Exception as e:
            return await state.end(
                warning=str(e).capitalize()
            )
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @commands.command(name="backup", help="DB backup.")
    @developer_predicator()
    async def text_backup(self, ctx: commands.Context):
        state = StateService(ctx)
        db = Database(directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except Exception as e:
            return await state.end(
                warning=str(e).capitalize()
            )
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    @app_commands.command(name="cogs", description="Lists cogs.")
    @developer_predicator()
    async def list_cogs_app_command(self, interaction: discord.Interaction):
        state = StateService(interaction)
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"Cogs for {interaction.guild.me.name}",
            color=discord.Color.blurple(),
        )
        for cog in sorted(DISCORD_COGS_CLASSES):
            if cog in self.bot.cogs:
                loaded.append(cog)
            else:
                not_loaded.append(cog)
        if loaded:
            embed.add_field(name="Loaded", value="\n".join(loaded), inline=False)
        if not_loaded:
            embed.add_field(
                name="Not Loaded", value="\n".join(not_loaded), inline=False
            )
        if not loaded and not not_loaded:
            embed.add_field(name="No cogs available.", inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    @commands.command(name="cogs", help="Lists cogs.")
    @developer_predicator()
    async def list_cogs_text_command(self, ctx: commands.Context):
        state = StateService(ctx)
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"Cogs for {ctx.guild.me.name}",
            color=discord.Color.blurple(),
        )
        for cog in sorted(DISCORD_COGS_CLASSES):
            if cog in self.bot.cogs:
                loaded.append(cog)
            else:
                not_loaded.append(cog)
        if loaded:
            embed.add_field(name="Loaded", value="\n".join(loaded), inline=False)
        if not_loaded:
            embed.add_field(
                name="Not Loaded", value="\n".join(not_loaded), inline=False
            )
        if not loaded and not not_loaded:
            embed.add_field(name="No cogs available.", inline=False)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    @app_commands.command(
        name="dlog", description="Resolve or update the notes on an issue by reference."
    )
    @app_commands.describe(
        reference="Specify the issue reference ID.",
        action="'resolve' or 'append' or 'overwrite'.",
        notes="Optionally specify notes to append or overwrite.",
    )
    @developer_predicator()
    async def update_developer_logs_app_command(
        self,
        ctx: commands.Context,
        reference: str,
        action: str,
        notes: Optional[str] = None,
    ):
        state = StateService(ctx)
        developer_log = await DeveloperLog.select(id=reference, resolved=False)
        if not developer_log:
            try:
                return await state.end(
                    warning=f"Issue not found. Received: {reference}."
                )
            except Exception as e:
                return await state.end(error=str(e).capitalize())
        if action and action.lower() == "resolve":
            await developer_log.resolve()
            detail = (
                "resolved the issue. The record "
                "will remain in the database for the next 30 days."
            )
        elif action and action.lower() == "append":
            await developer_log.append(notes)
            detail = "appended to the previous notes."
        elif action and action.lower() == "overwrite":
            await developer_log.overwrite(notes)
            detail = "overwrote the previous notes."
        try:
            return await state.end(
                success=f"\U000026a0\U0000fe0f " f"You successfully {detail}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    @commands.command(
        name="dlog", help="Resolve or update the notes on an issue by reference"
    )
    @developer_predicator()
    async def update_developer_logs_text_command(
        self,
        ctx: commands.Context,
        reference: str = commands.parameter(
            default=None, description="Specify the developer log reference ID."
        ),
        action: str = commands.parameter(
            default="append",
            description="Specify one of: `resolve` or `append` or `overwrite`.",
        ),
        *,
        notes: Optional[str] = commands.parameter(
            default=None, description="Optionally specify notes."
        ),
    ):
        state = StateService(ctx)
        developer_log = await DeveloperLog.select(id=reference, resolved=False)
        if not developer_log:
            try:
                return await state.end(
                    warning=f"Issue not found. Received: {reference}."
                )
            except Exception as e:
                return await state.end(error=str(e).capitalize())
        if action and action.lower() == "resolve":
            await developer_log.resolve()
            detail = "resolved the issue"
        elif action and action.lower() == "append":
            await developer_log.append(notes)
            detail = "appended to the previous notes"
        elif action and action.lower() == "overwrite":
            await developer_log.overwrite(notes)
            detail = "overwrote the previous notes"
        try:
            return await state.end(
                success=f"\U000026a0\U0000fe0f " f"You successfully {detail}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    @app_commands.command(name="dlogs", description="List issues.")
    @app_commands.describe(
        target="Specify one of: `all`, server ID or UUID.",
        filter="Optionally specify `resolved` or `unresolved`.",
    )
    @developer_predicator()
    async def list_developer_logs_app_command(
        self,
        interaction: discord.Interaction,
        target: str = None,
        filter: str = None,
    ):
        state = StateService(interaction)

        try:
            target_uuid = UUID(str(target))
            developer_logs = await DeveloperLog.select(id=target_uuid)
            title = f"{get_random_emoji()} Developer Logs"
        except (ValueError, TypeError) as e:
            logger.warning(str(e).capitalize())
            developer_logs, title = await resolve_objects(
                ctx_interaction_or_message=interaction,
                obj=DeveloperLog,
                state=state,
                target=target,
            )

        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            logger.warning(str(e).capitalize())
            is_at_home = False

        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}

        for developer_log in developer_logs:
            guild_dictionary.setdefault(developer_log.guild_snowflake, {"messages": {}})
            messages = guild_dictionary[developer_log.guild_snowflake]["messages"]
            messages.setdefault(
                developer_log.message_snowflake,
                {
                    "channel_snowflake": developer_log.channel_snowflake,
                    "developer_snowflakes": [],
                    "id": developer_log.id,
                    "notes": [],
                    "resolved": developer_log.resolved,
                },
            )
            messages[developer_log.message_snowflake]["developer_snowflakes"].extend(
                developer_log.developer_snowflakes
            )
            messages[developer_log.message_snowflake]["notes"].append(
                developer_log.notes
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_messages = await generate_skipped_messages(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_messages=skipped_messages,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for message_snowflake, entry in guild_data.get("messages", {}).items():
                channel = self.bot.get_channel(entry["channel_snowflake"])
                if not channel:
                    continue
                if filter == "resolved" and not entry.get("resolved"):
                    continue
                if filter == "unresolved" and entry.get("resolved"):
                    continue
                lines = []
                msg = await channel.fetch_message(message_snowflake)
                lines.append(
                    f'**Resolved:** {"\u2705" if entry.get("resolved") else "\u274c"}'
                )
                lines.append(f"**Message:** {msg.jump_url}")
                if target and str(target) == str(entry["id"]):
                    lines.append(
                        f'**Notes:** {entry["notes"] if entry.get("notes") is not None else None}'
                    )
                    lines.append(
                        f'**Assigned to:** {", ".join(str(d) for d in entry["developer_snowflakes"]) if entry.get("developer_snowflakes") else None}'
                    )
                else:
                    lines.append(f'**Reference:** {entry["id"]}')
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name=f"**Channel:** {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name=f"**Channel:** {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_messages:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_messages,
                    title="Skipped Messages in Server",
                )

        await StateService.send_pages(obj=DeveloperLog, pages=pages, state=state)

    @commands.command(name="dlogs", help="List issues.")
    @developer_predicator()
    async def list_developer_logs_text_command(
        self,
        ctx: commands.Context,
        *,
        target: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: `all`, server ID or UUID.",
        ),
        filter: Optional[str] = commands.parameter(
            default=None,
            description="Optionally specify `resolved` or `unresolved`.",
        ),
    ):
        state = StateService(ctx)

        try:
            target_uuid = UUID(str(target))
            developer_logs = await DeveloperLog.select(id=target_uuid)
            title = f"{get_random_emoji()} Developer Logs"
        except (ValueError, TypeError) as e:
            logger.warning(str(e).capitalize())
            developer_logs, title = await resolve_objects(
                ctx_interaction_or_message=ctx,
                obj=DeveloperLog,
                state=state,
                target=target,
            )

        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            logger.warning(str(e).capitalize())
            is_at_home = False

        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}

        for developer_log in developer_logs:
            guild_dictionary.setdefault(developer_log.guild_snowflake, {"messages": {}})
            messages = guild_dictionary[developer_log.guild_snowflake]["messages"]
            messages.setdefault(
                developer_log.message_snowflake,
                {
                    "channel_snowflake": developer_log.channel_snowflake,
                    "developer_snowflakes": [],
                    "id": developer_log.id,
                    "notes": [],
                    "resolved": developer_log.resolved,
                },
            )
            messages[developer_log.message_snowflake]["developer_snowflakes"].extend(
                developer_log.developer_snowflakes
            )
            messages[developer_log.message_snowflake]["notes"].append(
                developer_log.notes
            )

        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_messages = await generate_skipped_messages(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(
            guild_dictionary=guild_dictionary,
            skipped_guilds=skipped_guilds,
            skipped_messages=skipped_messages,
        )

        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for message_snowflake, entry in guild_data.get("messages", {}).items():
                channel = self.bot.get_channel(entry["channel_snowflake"])
                if not channel:
                    continue
                if filter == "resolved" and not entry.get("resolved"):
                    continue
                if filter == "unresolved" and entry.get("resolved"):
                    continue
                lines = []
                msg = await channel.fetch_message(message_snowflake)
                lines.append(
                    f'**Resolved:** {"\u2705" if entry.get("resolved") else "\u274c"}'
                )
                lines.append(f"**Message:** {msg.jump_url}")
                if target and str(target) == str(entry["id"]):
                    lines.append(
                        f'**Notes:** {entry["notes"] if entry.get("notes") is not None else None}'
                    )
                    lines.append(
                        f'**Assigned to:** {", ".join(str(d) for d in entry["developer_snowflakes"]) if entry.get("developer_snowflakes") else None}'
                    )
                else:
                    lines.append(f'**Reference:** {entry["id"]}')
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name=f"**Channel:** {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name=f"**Channel:** {channel.mention}",
                    value="\n".join(lines),
                    inline=False,
                )
            pages.append(embed)

        if is_at_home:
            if skipped_guilds:
                pages = generate_skipped_set_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_guilds,
                    title="Skipped Servers",
                )
            if skipped_messages:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_messages,
                    title="Skipped Messages in Server",
                )

        await StateService.send_pages(obj=DeveloperLog, pages=pages, state=state)

    # DONE
    @app_commands.command(
        name="load", description="Loads a cog by name 'vyrtuous.cog.<cog_name>.'"
    )
    @developer_predicator()
    async def load_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = StateService(interaction)
        try:
            await interaction.client.load_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(
                    warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=str(e).capitalize())
        try:
            return await state.end(
                success=f"Successfully loaded {module}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @commands.command(
        name="load", help="Loads a cog by name 'vyrtuous.cog.<cog_name>.'"
    )
    @developer_predicator()
    async def load_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(ctx)
        try:
            await self.bot.load_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(
                    warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=str(e).capitalize())
        try:
            return await state.end(
                success=f"Successfully loaded {module}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @app_commands.command(name="ping", description="Ping me!")
    @developer_predicator()
    async def ping_app_command(self, interaction: discord.Interaction):
        state = StateService(interaction)
        try:
            return await state.end(success=f"Pong!")
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @commands.command(name="ping", help="Ping me!")
    @developer_predicator()
    async def ping_text_command(self, ctx: commands.Context):
        state = StateService(ctx)
        try:
            return await state.end(success=f"Pong!")
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @app_commands.command(
        name="reload", description="Reloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @app_commands.check(at_home)
    @developer_predicator()
    async def reload_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = StateService(interaction)
        try:
            await interaction.client.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(
                    warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=str(e).capitalize())
        try:
            return await state.end(
                success=f"Successfully reloaded {module}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @commands.command(
        name="reload", help="Reloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @developer_predicator()
    async def reload_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(
                    warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=str(e).capitalize())
        try:
            return await state.end(
                success=f"Successfully reloaded {module}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @app_commands.command(name="sync", description="Sync app commands.")
    @developer_predicator()
    async def sync_app_command(
        self,
        interaction: discord.Interaction,
        spec: Optional[Literal["~", "*", "^"]] = None,
    ):
        await interaction.response.defer(ephemeral=True)
        state = StateService(interaction)
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
                try:
                    return await state.end(
                        warning=str(e).capitalize()
                    )
                except Exception as e:
                    return await state.end(error=str(e).capitalize())
        ret = 0
        for guild in guilds:
            try:
                await interaction.client.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        try:
            return await state.end(
                success=f"Synced the tree to {ret}/{len(guilds)}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @commands.command(name="sync", help="Sync app commands.")
    @developer_predicator()
    async def sync_text_command(
        self,
        ctx: commands.Context,
        guilds: commands.Greedy[discord.Object],
        spec: Optional[Literal["~", "*", "^"]] = None,
    ):
        state = StateService(ctx)
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
                try:
                    return await state.end(
                        warning=str(e).capitalize()
                    )
                except Exception as e:
                    return await state.end(error=str(e).capitalize())
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        try:
            return await state.end(
                success=f"Synced the tree to {ret}/{len(guilds)}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @app_commands.command(
        name="unload", description="Unloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @developer_predicator()
    async def unload_app_command(self, interaction: discord.Interaction, module: str):
        await interaction.response.defer(ephemeral=True)
        state = StateService(interaction)
        try:
            await interaction.client.unload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(
                    warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=str(e).capitalize())
        try:
            return await state.end(
                success=f"Successfully unloaded {module}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())

    # DONE
    @commands.command(
        name="unload", help="Unloads a cog by name 'vyrtuous.cog.<cog_name>'."
    )
    @developer_predicator()
    async def unload_text_command(self, ctx: commands.Context, *, module: str):
        state = StateService(ctx)
        try:
            await self.bot.reload_extension(module)
        except commands.ExtensionError as e:
            try:
                return await state.end(
                    warning=f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=str(e).capitalize())
        try:
            return await state.end(
                success=f"Successfully unloaded {module}."
            )
        except Exception as e:
            return await state.end(error=str(e).capitalize())


async def setup(bot: DiscordBot):
    await bot.add_cog(DevCommands(bot))
