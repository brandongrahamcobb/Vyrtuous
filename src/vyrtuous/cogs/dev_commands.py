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

from discord import app_commands
from discord.ext import commands
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import DISCORD_COGS_CLASSES
from vyrtuous.service.check_service import (
    at_home,
    developer_predicator
)
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.scope_service import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    resolve_objects,
    generate_skipped_channels,
    generate_skipped_guilds,
    generate_skipped_messages,
    clean_guild_dictionary,
    flush_page
)
from vyrtuous.service.state_service import StateService
from vyrtuous.database.database import Database
from vyrtuous.database.logs.developer_log import DeveloperLog
from vyrtuous.utils.emojis import get_random_emoji
import discord
from vyrtuous.utils.setup_logging import logger

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
                warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
            )
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
            )
        try:
            return await state.end(success=discord.File(db.file_name))
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    @app_commands.command(name="cogs", description="Lists cogs.")
    @developer_predicator()
    async def list_cogs_app_command(self, interaction: discord.Interaction):
        state = StateService(interaction)
        loaded, not_loaded = [], []
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"Cogs for {interaction.guild.me.name}",
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
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                    warning=f"\U000026a0\U0000fe0f "
                    f"Issue not found. Received: {reference}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c " f"{str(e).capitalize()}")
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
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                    warning=f"\U000026a0\U0000fe0f "
                    f"Issue not found. Received: {reference}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
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
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    @app_commands.command(name="dlogs", description="List issues.")
    @app_commands.describe(
        target="Specify one of: 'all', 'resolved' or 'unresolved'",
        value="Specify one of: channel ID/mention, reference ID and server ID.",
    )
    @developer_predicator()
    async def list_developer_logs_app_command(
        self,
        interaction: discord.Interaction,
        target: Optional[str],
        value: Optional[str],
    ):
        state = StateService(interaction)
        developer_logs, title = await resolve_objects(ctx_interaction_or_message=interaction, obj=DeveloperLog, state=state, target=target)
        
        try:
            is_at_home = at_home(ctx_interaction_or_message=interaction)
        except Exception as e:
            logger.warning(f"{str(e).capitalize()}")
            pass
    
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}

        for developer_log in developer_logs:
            guild_dictionary.setdefault(developer_log.guild_snowflake, {})
            guild_dictionary[developer_log.guild_snowflake].setdefault(
                developer_log.channel_snowflake, []
            )
            guild_dictionary[developer_log.guild_snowflake][
                developer_log.channel_snowflake
            ].append(
                {
                    "developer_snowflakes": developer_log.developer_snowflakes,
                    "id": developer_log.id,
                    "message_snowflake": developer_log.message_snowflake,
                    "notes": developer_log.notes,
                    "resolved": developer_log.resolved,
                }
            )

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_messages = await generate_skipped_messages(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_messages=skipped_messages)

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_logs in channels.items():
                channel = guild.get_channel(channel_snowflake)
                lines = []
                for member_data in channel_logs:
                    msg = await channel.fetch_message(
                        member_data["message_snowflake"]
                    )
                    lines.append(f"**Message:** {msg.jump_url}")
                    if member_data["resolved"] == False:
                        resolved = "\u274c"
                    elif member_data["resolved"] == True:
                        resolved = "\u2705"
                    lines.append(f"{resolved}**Reference:** {member_data['id']}")
                    if value:
                        lines.append(f"**Notes:** {member_data['notes']}")
                    if member_data["developer_snowflakes"]:
                        lines.append(
                            f"**Assigned to:** "
                            f"{', '.join(member_data['developer_snowflakes'])}"
                        )
                    field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name=f"**Channel:** {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
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
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
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
            description="Specify one of: `all`, `resolved`` or `unresolved`.",
        ),
        value: Optional[str] = commands.parameter(
            default=None,
            description="Specify one of: channel ID/mention, "
            "reference ID or server ID.",
        ),
    ):
        state = StateService(ctx)
        developer_logs, title = await resolve_objects(ctx_interaction_or_message=ctx, obj=DeveloperLog, state=state, target=target)
        
        try:
            is_at_home = at_home(ctx_interaction_or_message=ctx)
        except Exception as e:
            logger.warning(f"{str(e).capitalize()}")
            pass
    
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}

        for developer_log in developer_logs:
            guild_dictionary.setdefault(developer_log.guild_snowflake, {})
            guild_dictionary[developer_log.guild_snowflake].setdefault(
                developer_log.channel_snowflake, []
            )
            guild_dictionary[developer_log.guild_snowflake][
                developer_log.channel_snowflake
            ].append(
                {
                    "developer_snowflakes": developer_log.developer_snowflakes,
                    "id": developer_log.id,
                    "message_snowflake": developer_log.message_snowflake,
                    "notes": developer_log.notes,
                    "resolved": developer_log.resolved,
                }
            )

        skipped_channels = generate_skipped_channels(guild_dictionary)
        skipped_guilds = generate_skipped_guilds(guild_dictionary)
        skipped_messages = await generate_skipped_messages(guild_dictionary)
        guild_dictionary = clean_guild_dictionary(guild_dictionary=guild_dictionary, skipped_channels=skipped_channels, skipped_guilds=skipped_guilds, skipped_messages=skipped_messages)

        for guild_snowflake, channels in guild_dictionary.items():
            field_count = 0
            guild = self.bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_logs in channels.items():
                channel = guild.get_channel(channel_snowflake)
                lines = []
                for member_data in channel_logs:
                    msg = await channel.fetch_message(
                            member_data["message_snowflake"]
                        )
                    lines.append(f"**Message:** {msg.jump_url}")
                    if member_data["resolved"] == False:
                        resolved = "\u274c"
                    elif member_data["resolved"] == True:
                        resolved = "\u2705"
                    lines.append(f"{resolved}**Reference:** " f"{member_data['id']}")
                    if value:
                        lines.append(f"**Notes:** " f"{member_data['notes']}")
                    if member_data["developer_snowflakes"]:
                        lines.append(
                            f"**Assigned to:** "
                            f"{', '.join(member_data['developer_snowflakes'])}"
                        )
                    field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name=f"**Channel:** {channel.mention}",
                        value="\n".join(lines),
                        inline=False,
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
            if lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
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
            if skipped_channels:
                pages = generate_skipped_dict_pages(
                    chunk_size=chunk_size,
                    field_count=field_count,
                    pages=pages,
                    skipped=skipped_channels,
                    title="Skipped Channels in Server",
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
                    warning=f"\U000026a0\U0000fe0f "
                    f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} "
                f"Successfully loaded {module}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                    warning=f"\U000026a0\U0000fe0f "
                    f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} "
                f"Successfully loaded {module}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @app_commands.command(name="ping", description="Ping me!")
    @developer_predicator()
    async def ping_app_command(self, interaction: discord.Interaction):
        state = StateService(interaction)
        try:
            return await state.end(success=f"{get_random_emoji()} " f"Pong!")
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    @commands.command(name="ping", help="Ping me!")
    @developer_predicator()
    async def ping_text_command(self, ctx: commands.Context):
        state = StateService(ctx)
        try:
            return await state.end(success=f"{get_random_emoji()} " f"Pong!")
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                    warning=f"\U000026a0\U0000fe0f "
                    f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} "
                f"Successfully reloaded {module}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                    warning=f"\U000026a0\U0000fe0f "
                    f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} "
                f"Successfully reloaded {module}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                return await state.end(
                    success=f"{get_random_emoji()} " f"{msg}"
                )
            except Exception as e:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
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
                success=f"{get_random_emoji()} "
                f"Synced the tree to {ret}/{len(guilds)}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                return await state.end(
                    success=f"{get_random_emoji()} " f"{msg}"
                )
            except Exception as e:
                try:
                    return await state.end(
                        warning=f"\U000026a0\U0000fe0f " f"{str(e).capitalize()}"
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")
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
                success=f"{get_random_emoji()} "
                f"Synced the tree to {ret}/{len(guilds)}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                    warning=f"\U000026a0\U0000fe0f "
                    f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} "
                f"Successfully unloaded {module}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

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
                    warning=f"\U000026a0\U0000fe0f "
                    f"{e.__class__.__name__}: {str(e).capitalize()}"
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            return await state.end(
                success=f"{get_random_emoji()} "
                f"Successfully unloaded {module}."
            )
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")


async def setup(bot: DiscordBot):
    await bot.add_cog(DevCommands(bot))
