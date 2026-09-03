"""!/bin/python3
coordinator_service.py The purpose of this program is to extend Service to service the coordinator class.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.permissions import PermissionGroup, PermissionScope
from vyrtuous.cache.registry import MemberState, PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.permission_entry import PermissionEntry
from vyrtuous.listing import list_service
from vyrtuous.permissions import permission_service
from vyrtuous.utils.errors.error import CheckFailure
from vyrtuous.utils.messaging import emojis

MODEL = PermissionEntry


async def build_dictionary(
    group: PermissionGroup,
    guild_snowflake: int | None,
    channel_snowflake: int | None,
) -> dict[int | str, dict[str, dict[int, dict[str, dict[int, dict[str, str]]]]]]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    records = []
    dictionary: dict[
        int | str, dict[str, dict[int, dict[str, dict[int, dict[str, str]]]]]
    ] = {}
    if group.scope is PermissionScope.GLOBAL:
        if guild_snowflake or channel_snowflake:
            return dictionary
        select_kwargs = {}
    elif group.scope is PermissionScope.GUILD:
        if channel_snowflake:
            return dictionary
        select_kwargs = {"guild_snowflake": guild_snowflake} if guild_snowflake else {}
    elif group.scope is PermissionScope.CHANNEL:
        select_kwargs = {}
        if guild_snowflake:
            select_kwargs["guild_snowflake"] = guild_snowflake
        if channel_snowflake:
            select_kwargs["channel_snowflake"] = channel_snowflake
    else:
        return dictionary
    records = await database_factory.select(
        group_alias=group.alias,
        singular=False,
        **select_kwargs,
    )
    for record in records:
        key: int | str = (
            "global"
            if group.scope is PermissionScope.GLOBAL
            else record.guild_snowflake
        )
        dictionary.setdefault(key, {"members": {}})
        members = dictionary[key]["members"]
        members.setdefault(record.member_snowflake, {"groups": {}})
        members[record.member_snowflake]["groups"].setdefault(
            record.channel_snowflake, {}
        )
        members[record.member_snowflake]["groups"][record.channel_snowflake].update(
            {"group_alias": group.alias}
        )
    return dictionary


async def build_pages(
    group: PermissionGroup,
    member_snowflake: int | None,
    guild_snowflake: int | None,
    channel_snowflake: int | None,
) -> str | list[discord.Embed]:
    bot: DiscordBot = DiscordBot.get_instance()
    lines: list[str] = []
    pages: list[discord.Embed] = []

    if member_snowflake:
        member = bot.get_user(member_snowflake)
        if member is None:
            simplified_member = bot.registry.get(MemberState).active.get(
                member_snowflake, None
            )
            if simplified_member:
                obj_name = simplified_member[0]
            else:
                return "This command must target a valid member."
        else:
            obj_name = member.display_name
    elif guild_snowflake and channel_snowflake:
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            return f"No {group.name}s found."
        channel = guild.get_channel(channel_snowflake)
        if channel is None:
            return f"No {group.name}s found."
        obj_name = channel.name
    elif guild_snowflake:
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            return f"No {group.name}s found."
        obj_name = guild.name
    else:
        return f"No {group.name}s found."
    title = f"{emojis.get_random_emoji()} {group.name}s for {obj_name}"
    dictionary = await build_dictionary(
        group=group,
        guild_snowflake=guild_snowflake,
        channel_snowflake=channel_snowflake,
    )
    for scope_key, scope_data in dictionary.items():
        member_n = 0
        field_count = 0
        lines: list[str] = []
        if isinstance(scope_key, str):
            guild = None
            description = "Global"
        else:
            guild = bot.get_guild(scope_key)
            if guild is None:
                continue
            description = guild.name
        embed = discord.Embed(
            title=title, description=description, color=discord.Color.blue()
        )
        for member_snowflake_, member_data in scope_data.get("members", {}).items():
            if member_snowflake is not None and member_snowflake_ != member_snowflake:
                continue
            member = (
                guild.get_member(member_snowflake_)
                if guild
                else bot.get_user(member_snowflake_)
            )
            if member:
                display = f"{member.mention}"
            else:
                simplified_member = bot.registry.get(MemberState).active.get(
                    member_snowflake_, None
                )
                if not simplified_member:
                    continue
                display = f"{simplified_member[0]} ({member_snowflake_})"
            member_n += 1
            groups = member_data.get("groups", {})
            if group.scope is not PermissionScope.CHANNEL:
                lines.append(f"**User:** {display}")
                field_count += 1
            else:
                for channel_snowflake_ in groups:
                    channel = guild.get_channel(channel_snowflake_) if guild else None
                    if channel is None:
                        continue
                    lines.append(f"**User:** {display}")
                    field_count += 1
            if field_count >= list_service.CHUNK_SIZE:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
                embed = list_service.flush_page(embed, pages, title, description)
                lines = []
                field_count = 0
        if lines:
            embed.add_field(name="Information", value="\n".join(lines), inline=False)
        if member_n == 0:
            continue
        original_description = embed.description or ""
        embed.description = f"**{original_description} ({member_n})**"
        pages.append(embed)
    if not pages:
        return f"No {group.name}s found."
    return pages

async def build_summary_pages(author_snowflake: int, display_name: str, guild_snowflake: int, member_snowflake: int):
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    all_assigned_groups = await permission_service.resolve_all_assigned_groups(permission_state=permission_state, member_snowflake=member_snowflake)
    lines = []
    global_groups = []
    guild_groups = {}
    channel_groups = {}
    for group, group_guild_snowflake, group_channel_snowflake in all_assigned_groups: 
        match group.scope:
            case PermissionScope.GLOBAL:
                try:
                     await permission_service.has_permissions_at_all(
                        permission_state=permission_state,
                        member_snowflake=author_snowflake,
                        requested=["command.info.scope.global"],
                    )
                except CheckFailure:
                    continue
                global_groups.append(group)
            case PermissionScope.GUILD:
                if group_guild_snowflake:
                    try:
                        await permission_service.has_permissions_at_all(
                            permission_state=permission_state,
                            member_snowflake=author_snowflake,
                            requested=["command.info.scope.guild"],
                        )
                    except CheckFailure:
                        continue
                    if guild_snowflake != group_guild_snowflake:
                        try:
                            await permission_service.has_permissions_at_all(
                                permission_state=permission_state,
                                member_snowflake=author_snowflake,
                                requested=["other_guilds"],
                            )
                        except CheckFailure:
                            continue
                else:
                    bot.logger.debug("Unexpected group with scope 'GUILD' for /groups command.")
                    return
                guild = bot.get_guild(group_guild_snowflake)
                if guild:
                    guild_groups.setdefault(guild.name, []).append(group)
                else:
                    continue
            case PermissionScope.CHANNEL:
                if group_channel_snowflake is not None and group_guild_snowflake is not None:
                    guild = bot.get_guild(group_guild_snowflake)
                    if guild is None:
                        continue
                    channel = guild.get_channel(group_channel_snowflake)
                    if channel is None:
                        continue
                    channel_groups.setdefault((guild.name, channel.name), []).append(group)
                else:
                    bot.logger.debug("Unexpected group with scope 'CHANNEL' for /groups command.")
                    return
    group_entries = []
    for group in global_groups:
        group_entries.append(('🌐 Global', None, f'• {group.name}'))
    for guild_name, groups in guild_groups.items():
        for group in groups:
            group_entries.append(('🏠 Server', guild_name, f'• {group.name}'))
    for (guild_name, channel_name), groups in channel_groups.items():
        for group in groups:
            group_entries.append(('💬 Channel', f'{guild_name} → #{channel_name}', f'• {group.name}'))
    pages = []
    embed = discord.Embed(color=discord.Color.blurple(), title=f'Groups for {display_name}')
    field_count = 0
    sections = {'🌐 Global': [], '🏠 Server': [], '💬 Channel': []}
    for section, parent, line in group_entries:
        sections[section].append((parent, line))
        field_count += 1
        if field_count >= list_service.CHUNK_SIZE:
            for section_name, entries in sections.items():
                if entries:
                    grouped = {}
                    for parent, entry in entries:
                        grouped.setdefault(parent, []).append(entry)
                    value = '\n'.join(
                        f'**{parent}**\n' + '\n'.join(values)
                        for parent, values in grouped.items()
                        if parent is not None
                    )
                    if section_name == '🌐 Global':
                        value = '\n'.join(values for _, values in entries)
                    embed.add_field(name=section_name, value=value or 'N/A', inline=False)
                else:
                    embed.add_field(name=section_name, value='N/A', inline=False)
            embed = list_service.flush_page(embed, pages, f'Groups for {display_name}', display_name)
            embed = discord.Embed(color=discord.Color.blurple(), title=f'Groups for {display_name}')
            sections = {'🌐 Global': [], '🏠 Server': [], '💬 Channel': []}
            field_count = 0
    if field_count:
        for section_name, entries in sections.items():
            if entries:
                grouped = {}
                for parent, entry in entries:
                    grouped.setdefault(parent, []).append(entry)
                value = '\n'.join(
                    f'**{parent}**\n' + '\n'.join(values)
                    for parent, values in grouped.items()
                    if parent is not None
                )
                if section_name == '🌐 Global':
                    value = '\n'.join(values for _, values in entries)
                embed.add_field(name=section_name, value=value or 'N/A', inline=False)
            else:
                embed.add_field(name=section_name, value='N/A', inline=False)
        pages.append(embed)
    return pages

