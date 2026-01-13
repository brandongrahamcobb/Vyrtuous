"""setup_logging.py The purpose of this program is to set up logging for Vyrtuous.

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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.channel_service import resolve_channel
from vyrtuous.service.check_service import role_check_without_specifics
from vyrtuous.service.member_service import resolve_member
from vyrtuous.service.role_service import resolve_role
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
    for guild_snowflake, list in skipped.items():
        embed = discord.Embed(
            color=discord.Color.red(), title=f"{title} ({guild_snowflake})"
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

async def resolve_where_kwargs(channel_obj, guild_obj, member_obj, role_obj):
    kwargs = {}
    match (channel_obj, guild_obj, member_obj, role_obj):
        case (c, g, None, None) if c and g:
            kwargs['channel_snowflake'] = c.id
            kwargs['guild_snowflake'] = g.id
        case (None, g, None, None) if g:
            kwargs['guild_snowflake'] = g.id
        case (None, g, m, None) if g and m:
            kwargs['guild_snowflake'] = g.id
            kwargs['member_snowflake'] = m.id
        case (None, g, None, r) if g and r:
            kwargs['guild_snowflake'] = g.id
            kwargs['role_snowflake'] = r.id
        case _:
            raise ValueError
    return kwargs

async def resolve_objects(
    ctx_interaction_or_message, obj, state, target
):
    channel_obj, guild_obj, member_obj, role_obj = None, None, None, None
    
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
            try:
                member_obj = await resolve_member(
                    ctx_interaction_or_message=ctx_interaction_or_message,
                    member_str=target,
                )
            except Exception as e:
                bot = DiscordBot.get_instance()
                guild_obj = bot.get_guild(int(target))
                if not guild_obj:
                    try:
                        role_obj = await resolve_role(
                            ctx_interaction_or_message=ctx_interaction_or_message,
                            role_str=target,
                        )
                    except Exception as e:
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
        role_obj=role_obj
    )
    objects = await obj.select(**where_kwargs)
    if not objects:
        try:
            msg = f"No {obj.PLURAL} exist for target: {target}."
            return await state.end(warning=f"\U000026a0\U0000fe0f {msg}")
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")
    
    title = f"{get_random_emoji()} {obj.PLURAL}"
    return objects, title
