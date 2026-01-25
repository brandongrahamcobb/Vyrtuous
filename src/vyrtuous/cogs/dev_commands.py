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
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.actions.hide import Hide
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.database import Database
from vyrtuous.db.mgmt.bug import Bug
from vyrtuous.db.roles.developer import developer_predicator
from vyrtuous.inc.helpers import DISCORD_COGS, DISCORD_COGS_CLASSES
from vyrtuous.fields.snowflake import (
    AppRoleSnowflake,
    RoleSnowflake,
)
from vyrtuous.utils.home import at_home
from vyrtuous.utils.logger import logger
from vyrtuous.service.message_service import MessageService
from vyrtuous.service.state_service import StateService
from vyrtuous.service.discord_object_service import (
    DiscordObject,
    DiscordObjectNotFound,
)
from vyrtuous.utils.guild_dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_guilds,
    generate_skipped_messages,
    clean_guild_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji


class DevCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @app_commands.command(
        name="assoc", description="Associate a hide or text-mute alias to a role."
    )
    @developer_predicator()
    @app_commands.describe(
        alias_name="One of: `hide`, `tmute`",
    )
    async def associate_alias_to_role_app_command(
        self, interaction: discord.Interaction, alias_name: str, role: AppRoleSnowflake
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        try:
            role_dict = await do.determine_from_target(target=role)
        except (DiscordObjectNotFound, TypeError) as e:
            logger.info(e)
            role = None
        else:
            role = interaction.guild.get_role(role_dict.get("id", None))
        alias = await Alias.select(
            alias_name=alias_name, guild_snowflake=interaction.guild.id, singular=True
        )
        if not alias:
            return await state.end(success=f"No such alias ({alias_name}) exists.")
        if alias and not getattr(alias, "role_snowflake"):
            if not role:
                for role in interaction.guild.roles:
                    if role.name == alias_name:
                        return await state.end(
                            warning=f"{role.name} already exists. You must specify it to override."
                        )
                role = await interaction.guild.create_role(name=alias_name)
                if alias.alias_type not in ("hide", "tmute"):
                    return await state.end(
                        warning=f"Alias `{alias.alias_name}` of type `{alias.alias_type}` "
                        f"cannot be associated with role {role.mention}."
                    )
        elif getattr(alias, "role_snowflake"):
            return await state.end(
                warning=f"Alias ({alias_name}) is already associated with a role."
            )
        channel = interaction.guild.get_channel(alias.channel_snowflake)
        try:
            if alias.alias_type == "hide":
                overwrite = discord.PermissionOverwrite(connect=False)
                try:
                    await channel.set_permissions(role, overwrite=overwrite)
                except discord.Forbidden as e:
                    logger.warning(e)
                hides = await Hide.select(
                    channel_snowflake=channel.id, guild_snowflake=interaction.guild.id
                )
                for hide in hides:
                    channel = interaction.guild.get_channel(hide.channel_snowflake)
                    member = interaction.guild.get_member(hide.member_snowflake)
                    try:
                        await channel.set_permissions(member, overwrite=None)
                    except discord.Forbidden as e:
                        logger.warning(e)
                    await Hide.administer_role(
                        guild_snowflake=interaction.guild.id,
                        member_snowflake=member.id,
                        role_snowflake=role.id,
                        state=state,
                    )
            elif alias.alias_type == "tmute":
                overwrite = discord.PermissionOverwrite(
                    send_messages=False,
                    add_reactions=False,
                    create_public_threads=False,
                    create_private_threads=False,
                    send_messages_in_threads=False,
                )
                try:
                    await channel.set_permissions(role, overwrite=overwrite)
                except discord.Forbidden as e:
                    logger.warning(e)
                text_mutes = await TextMute.select(
                    channel_snowflake=channel.id, guild_snowflake=interaction.guild.id
                )
                for text_mute in text_mutes:
                    channel = interaction.guild.get_channel(text_mute.channel_snowflake)
                    member = interaction.guild.get_member(text_mute.member_snowflake)
                    try:
                        await channel.set_permissions(member, overwrite=None)
                        await member.add_roles(
                            role, reason="Associating old text mutes with a role"
                        )
                    except discord.Forbidden as e:
                        logger.warning(e)
        except discord.Forbidden as e:
            logger.warning(e)
        where_kwargs = {
            "alias_name": alias.alias_name,
            "guild_snowflake": interaction.guild.id,
        }
        set_kwargs = {"role_snowflake": role.id}
        await Alias.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)

        return await state.end(
            success=f"Alias `{alias.alias_name}` of type `{alias.alias_type}` "
            f"was associated with role {role.mention}."
        )

    @commands.command(
        name="assoc", help="Associate a hide or text-mute alias to a role."
    )
    @developer_predicator()
    async def associate_alias_to_role_text_command(
        self,
        ctx: commands.Context,
        alias_name: str = commands.parameter(description="Alias/Pseudonym"),
        role: RoleSnowflake = commands.parameter(
            default=None, description="Tag a role or include the ID."
        ),
    ):

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        try:
            role_dict = await do.determine_from_target(target=role)
        except (DiscordObjectNotFound, TypeError) as e:
            logger.info(e)
            role = None
        else:
            role = ctx.guild.get_role(role_dict.get("id", None))
        alias = await Alias.select(
            alias_name=alias_name, guild_snowflake=ctx.guild.id, singular=True
        )
        if not alias:
            return await state.end(warning=f"No such alias ({alias_name}) exists.")
        if alias and not getattr(alias, "role_snowflake"):
            if not role:
                for role in ctx.guild.roles:
                    if role.name == alias_name:
                        return await state.end(
                            warning=f"{role.name} already exists. You must specify it to override."
                        )
                try:
                    role = await ctx.guild.create_role(name=alias_name)
                except discord.Forbidden as e:
                    logger.warning(e)
                if alias.alias_type not in ("hide", "tmute"):
                    return await state.end(
                        warning=f"Alias `{alias.alias_name}` of type `{alias.alias_type}` "
                        f"cannot be associated with role {role.mention}."
                    )
        elif getattr(alias, "role_snowflake"):
            return await state.end(
                warning=f"Alias ({alias_name}) is already associated with a role."
            )
        channel = ctx.guild.get_channel(alias.channel_snowflake)
        try:
            if alias.alias_type == "hide":
                overwrite = discord.PermissionOverwrite(connect=False)
                try:
                    await channel.set_permissions(role, overwrite=overwrite)
                except discord.Forbidden as e:
                    logger.warning(e)
                hides = await Hide.select(
                    channel_snowflake=channel.id, guild_snowflake=ctx.guild.id
                )
                for hide in hides:
                    channel = ctx.guild.get_channel(hide.channel_snowflake)
                    member = ctx.guild.get_member(hide.member_snowflake)
                    try:
                        await channel.set_permissions(member, overwrite=None)
                    except discord.Forbidden as e:
                        logger.warning(e)
                    await Hide.administer_role(
                        guild_snowflake=ctx.guild.id,
                        member_snowflake=member.id,
                        role_snowflake=role.id,
                        state=state,
                    )
            elif alias.alias_type == "tmute":
                overwrite = discord.PermissionOverwrite(
                    send_messages=False,
                    add_reactions=False,
                    create_public_threads=False,
                    create_private_threads=False,
                    send_messages_in_threads=False,
                )
                try:
                    await channel.set_permissions(role, overwrite=overwrite)
                except discord.Forbidden as e:
                    logger.warning(e)
                text_mutes = await TextMute.select(
                    channel_snowflake=channel.id, guild_snowflake=ctx.guild.id
                )
                for text_mute in text_mutes:
                    channel = ctx.guild.get_channel(text_mute.channel_snowflake)
                    member = ctx.guild.get_member(text_mute.member_snowflake)
                    try:
                        await channel.set_permissions(member, overwrite=None)
                        await member.add_roles(
                            role, reason="Associating old text mutes with the new "
                        )
                    except discord.Forbidden as e:
                        logger.warning(e)
        except discord.Forbidden as e:
            logger.warning(e)

        where_kwargs = {"alias_name": alias.alias_name, "guild_snowflake": ctx.guild.id}
        set_kwargs = {"role_snowflake": role.id}
        await Alias.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)

        return await state.end(
            success=f"Alias `{alias.alias_name}` of type `{alias.alias_type}` "
            f"was associated with role {role.mention}."
        )

    # DONE
    @app_commands.command(name="backup", description="DB backup.")
    @developer_predicator()
    async def app_backup(self, interaction: discord.Interaction):
        await interaction.response.defer(ephemeral=True)
        state = StateService(source=interaction)
        db = Database(directory="/app/backups")
        try:
            db.create_backup_directory()
            db.execute_backup()
        except RuntimeError as e:
            return await state.end(warning=str(e).capitalize())
        return await state.end(success=discord.File(db.file_name))

    # DONE
    @commands.command(name="backup", help="DB backup.")
    @developer_predicator()
    async def text_backup(self, ctx: commands.Context):
        state = StateService(source=ctx)
        db = Database(directory="/app/backups")
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
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=interaction)
        title = f"{get_random_emoji()} Developer Logs"

        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        try:
            target_uuid = UUID(str(target))
            kwargs = {"id": target_uuid}
        except Exception as e:
            logger.warning(str(e).capitalize())
            object_dict = await do.determine_from_target(target=target)
            kwargs = object_dict.get("columns", None)

        bugs = await Bug.select(**kwargs)

        for bug in bugs:
            guild_dictionary.setdefault(bug.guild_snowflake, {"messages": {}})
            messages = guild_dictionary[bug.guild_snowflake]["messages"]
            messages.setdefault(
                bug.message_snowflake,
                {
                    "channel_snowflake": bug.channel_snowflake,
                    "developer_snowflakes": [],
                    "id": bug.id,
                    "notes": [],
                    "resolved": bug.resolved,
                },
            )
            messages[bug.message_snowflake]["developer_snowflakes"].extend(
                bug.member_snowflakes
            )
            messages[bug.message_snowflake]["notes"].append(bug.notes)

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

        await StateService.send_pages(obj=Bug, pages=pages, state=state)

    @commands.command(name="dlogs", help="List issues.")
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
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        is_at_home = at_home(source=ctx)
        title = f"{get_random_emoji()} Developer Logs"

        state = StateService(source=ctx)
        do = DiscordObject(ctx=ctx)
        try:
            target_uuid = UUID(str(target))
            kwargs = {"id": target_uuid}
        except Exception as e:
            logger.warning(str(e).capitalize())
            object_dict = await do.determine_from_target(target=target)
            kwargs = object_dict.get("columns", None)

        bugs = await Bug.select(**kwargs)

        for bug in bugs:
            guild_dictionary.setdefault(bug.guild_snowflake, {"messages": {}})
            messages = guild_dictionary[bug.guild_snowflake]["messages"]
            messages.setdefault(
                bug.message_snowflake,
                {
                    "channel_snowflake": bug.channel_snowflake,
                    "developer_snowflakes": [],
                    "id": bug.id,
                    "notes": [],
                    "resolved": bug.resolved,
                },
            )
            messages[bug.message_snowflake]["developer_snowflakes"].extend(
                bug.member_snowflakes
            )
            messages[bug.message_snowflake]["notes"].append(bug.notes)

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

        await StateService.send_pages(obj=Bug, pages=pages, state=state)

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
        guilds: commands.Greedy[discord.Object],
        spec: Optional[Literal["~", "*", "^"]] = None,
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
