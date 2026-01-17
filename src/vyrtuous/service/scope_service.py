"""scope_service.py The purpose of this program is to manage list command logic.

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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.check_service import role_check_without_specifics
from vyrtuous.service.resolution.channel_service import resolve_channel
from vyrtuous.service.resolution.member_service import resolve_member
from vyrtuous.service.resolution.role_service import resolve_role
from vyrtuous.service.logging_service import logger
from vyrtuous.utils.emojis import get_random_emoji


def generate_skipped_set_pages(chunk_size, field_count, pages, skipped, title):
    embed = discord.Embed(title=title, description="\u200b", color=discord.Color.blue())
    lines = []
    for snowflake in skipped:
        if field_count >= chunk_size:
            embed.description = "\n".join(lines)
            pages.append(embed)
            embed = discord.Embed(
                title=f"{title} continued...", color=discord.Color.red()
            )
            lines = []
            field_count = 0
        lines.append(str(snowflake))
        field_count += 1
    embed.description = "\n".join(lines)
    pages.append(embed)
    return pages


def generate_skipped_dict_pages(chunk_size, field_count, pages, skipped, title):
    bot = DiscordBot.get_instance()
    for guild_snowflake, list in skipped.items():
        guild = bot.get_guild(guild_snowflake)
        embed = discord.Embed(
            color=discord.Color.red(), title=f"{title} ({guild.name})"
        )
        field_count = 0
        lines = []
        for snowflake in list:
            if field_count >= chunk_size:
                embed.description = "\n".join(lines)
                pages.append(embed)
                embed = discord.Embed(
                    color=discord.Color.red(),
                    title=f"{title} ({guild_snowflake}) continued...",
                )
                field_count = 0
                lines = []
            lines.append(str(snowflake))
            field_count += 1
        embed.description = "\n".join(lines)
        pages.append(embed)
    return pages


async def resolve_where_kwargs(channel_obj, guild_obj, member_obj, role_obj):
    kwargs = {}
    match (channel_obj, guild_obj, member_obj, role_obj):
        case (c, None, None, None) if c:
            kwargs["channel_snowflake"] = c.id
        case (None, g, None, None) if g:
            kwargs["guild_snowflake"] = g.id
        case (None, None, m, None) if m:
            kwargs["member_snowflake"] = m.id
        case (None, None, None, r) if r:
            kwargs["role_snowflake"] = r.id
        case (c, g, None, None) if c and g:
            kwargs["channel_snowflake"] = c.id
            kwargs["guild_snowflake"] = g.id
        case (None, g, None, None) if g:
            kwargs["guild_snowflake"] = g.id
        case (c, g, m, None) if c and g and m:
            kwargs["channel_snowflake"] = c.id
            kwargs["guild_snowflake"] = g.id
            kwargs["member_snowflake"] = m.id
        case (None, g, m, None) if g and m:
            kwargs["guild_snowflake"] = g.id
            kwargs["member_snowflake"] = m.id
        case (c, None, None, None) if c:
            kwargs["channel_snowflake"] = c.id
        case (None, g, None, r) if g and r:
            kwargs["guild_snowflake"] = g.id
            kwargs["role_snowflake"] = r.id
        case (None, None, None, None):
            return {}
        case _:
            raise ValueError
    return kwargs


async def resolve_objects(ctx_interaction_or_message, obj, state, target):
    channel_obj, guild_obj, member_obj, role_obj = None, None, None, None
    title = f"{get_random_emoji()} {obj.PLURAL}"
    highest_role = await role_check_without_specifics(
        ctx_interaction_or_message=ctx_interaction_or_message
    )
    if target and target.lower() == "all":
        if highest_role not in ("System Owner", "Developer"):
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f "
                    f"You are not authorized to list {obj.PLURAL.lower()} "
                    f"across all servers."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
    elif target:
        try:
            channel_obj = await resolve_channel(
                ctx_interaction_or_message=ctx_interaction_or_message,
                channel_str=target,
            )
        except Exception as e:
            logger.info(str(e).capitalize())
            try:
                member_obj = await resolve_member(
                    ctx_interaction_or_message=ctx_interaction_or_message,
                    member_str=target,
                )
            except Exception as e:
                logger.info(str(e).capitalize())
                bot = DiscordBot.get_instance()
                guild_obj = bot.get_guild(int(target))
                if not guild_obj:
                    try:
                        role_obj = await resolve_role(
                            ctx_interaction_or_message=ctx_interaction_or_message,
                            role_str=target,
                        )
                    except Exception as e:
                        logger.info(str(e).capitalize())
                        try:
                            return await state.end(
                                warning=f"\U000026a0\U0000fe0f "
                                f"Scope must be one of: `all`, channel ID/mention, "
                                f"member ID/mention, role ID/mention, server ID or empty. Received: {target}."
                            )
                        except Exception as e:
                            return await state.end(
                                error=f"\u274c {str(e).capitalize()}"
                            )
                else:
                    if highest_role in ("Coordinator", "Moderator", "Everyone"):
                        try:
                            return await state.end(
                                warning=f"\U000026a0\U0000fe0f "
                                f"You are not authorized to list {obj.PLURAL.lower()} "
                                f"for specific servers."
                            )
                        except Exception as e:
                            return await state.end(
                                error=f"\u274c {str(e).capitalize()}"
                            )
    where_kwargs = await resolve_where_kwargs(
        channel_obj=channel_obj,
        guild_obj=guild_obj,
        member_obj=member_obj,
        role_obj=role_obj,
    )
    objects = await obj.select(**where_kwargs)
    return objects, title


def generate_skipped_guilds(guild_dictionary: dict) -> set:
    bot = DiscordBot.get_instance()
    skipped_guilds = set()
    for guild_snowflake in guild_dictionary:
        if not bot.get_guild(guild_snowflake):
            skipped_guilds.add(guild_snowflake)
    return skipped_guilds


def generate_skipped_channels(guild_dictionary: dict) -> dict:
    bot = DiscordBot.get_instance()
    skipped_channels = {}
    for guild_snowflake, guild_data in guild_dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if not guild:
            continue
        for channel_snowflake in guild_data.get("channels", {}):
            if not guild.get_channel(channel_snowflake):
                skipped_channels.setdefault(guild_snowflake, []).append(
                    channel_snowflake
                )
    return skipped_channels


def generate_skipped_members(guild_dictionary: dict) -> dict:
    bot = DiscordBot.get_instance()
    skipped_members = {}
    for guild_snowflake, guild_data in guild_dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if not guild:
            continue
        for member_snowflake in guild_data.get("members", {}):
            if not guild.get_member(member_snowflake):
                skipped_members.setdefault(guild_snowflake, []).append(member_snowflake)
    return skipped_members


def generate_skipped_roles(guild_dictionary: dict) -> dict:
    bot = DiscordBot.get_instance()
    skipped_roles = {}
    for guild_snowflake, guild_data in guild_dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if not guild:
            continue
        for role_snowflake in guild_data.get("roles", {}):
            if not guild.get_role(role_snowflake):
                skipped_roles.setdefault(guild_snowflake, []).append(role_snowflake)
    return skipped_roles


def generate_skipped_snowflakes(guild_dictionary: dict) -> dict:
    bot = DiscordBot.get_instance()
    skipped_snowflakes = {}
    for guild_snowflake, guild_data in guild_dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if not guild:
            continue
        for snowflake in guild_data.get("snowflakes", []):
            if not guild.get_channel(snowflake):
                if not guild.get_member(snowflake):
                    if not guild.get_role(snowflake):
                        skipped_snowflakes.append(snowflake)
    return skipped_snowflakes


def clean_guild_dictionary(
    guild_dictionary: dict,
    *,
    skipped_guilds: set | None = None,
    skipped_channels: dict | None = None,
    skipped_members: dict | None = None,
    skipped_messages: dict | None = None,
    skipped_roles: dict | None = None,
    skipped_snowflakes: dict | None = None,
) -> dict:
    cleaned = {}
    skipped_guilds = skipped_guilds or set()
    skipped_channels = skipped_channels or {}
    skipped_members = skipped_members or {}
    skipped_messages = skipped_messages or {}
    skipped_roles = skipped_roles or {}
    skipped_snowflakes = skipped_snowflakes or {}
    for guild_snowflake, guild_data in guild_dictionary.items():
        if guild_snowflake in skipped_guilds:
            continue
        channels = {
            channel_snowflake: entries
            for channel_snowflake, entries in guild_data.get("channels", {}).items()
            if channel_snowflake not in skipped_channels.get(guild_snowflake, [])
        }
        members = {
            member_snowflake: entries
            for member_snowflake, entries in guild_data.get("members", {}).items()
            if member_snowflake not in skipped_members.get(guild_snowflake, [])
        }
        messages = {
            message_snowflake: entries
            for message_snowflake, entries in guild_data.get("messages", {}).items()
            if message_snowflake not in skipped_messages.get(guild_snowflake, [])
        }
        roles = {
            role_snowflake: entries
            for role_snowflake, entries in guild_data.get("roles", {}).items()
            if role_snowflake not in skipped_roles.get(guild_snowflake, [])
        }
        snowflakes = {
            snowflake: entries
            for snowflake, entries in guild_data.get("snowflakes", {}).items()
            if snowflake not in skipped_snowflakes.get(guild_snowflake, [])
        }
        cleaned[guild_snowflake] = {
            "channels": channels,
            "members": members,
            "messages": messages,
            "roles": roles,
            "snowflakes": snowflakes,
        }
    return cleaned


async def generate_skipped_messages(guild_dictionary: dict) -> dict:
    bot = DiscordBot.get_instance()
    skipped_messages = {}
    for guild_snowflake, guild_data in guild_dictionary.items():
        guild = bot.get_guild(guild_snowflake)
        if not guild:
            continue
        for channel_snowflake, channel_logs in guild_data.get("channels", {}).items():
            channel = guild.get_channel(channel_snowflake)
            if not channel:
                continue
            for member_data in channel_logs:
                try:
                    await channel.fetch_message(member_data["message_snowflake"])
                except Exception:
                    skipped_messages.setdefault(guild_snowflake, []).append(
                        member_data["message_snowflake"]
                    )
    return skipped_messages


def flush_page(embed, pages, title, guild_name):
    pages.append(embed)
    return (
        discord.Embed(
            title=title,
            description=f"{guild_name} continued...",
            color=discord.Color.blue(),
        ),
        0,
    )
