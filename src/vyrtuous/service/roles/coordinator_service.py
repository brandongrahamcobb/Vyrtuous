"""coordinator.py The purpose of this program is to inherit from the DatabaseFactory to provide the coordinator role.
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

from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.roles.coordinator import Coordinator
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.service.roles.administrator_service import is_administrator_wrapper
from vyrtuous.service.roles.developer_service import is_developer_wrapper
from vyrtuous.service.roles.guild_owner_service import is_guild_owner_wrapper
from vyrtuous.service.roles.sysadmin_service import is_sysadmin_wrapper
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    generate_skipped_set_pages,
)
from vyrtuous.utils.dir_to_classes import skip_db_discovery
from vyrtuous.utils.emojis import get_random_emoji


@skip_db_discovery
class NotCoordinator(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a sysadmin, developer, guild owner, administrator or coordinator in this channel.",
    ):
        super().__init__(message)


async def is_coordinator_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_coordinator(
        channel_snowflake=source.channel.id,
        guild_snowflake=source.guild.id,
        member_snowflake=int(member_snowflake),
    )


async def is_coordinator(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    coordinator = await Coordinator.select(
        channel_snowflake=int(channel_snowflake),
        guild_snowflake=int(guild_snowflake),
        member_snowflake=int(member_snowflake),
        singular=True,
    )
    if not coordinator:
        raise NotCoordinator
    return True


async def is_coordinator_at_all(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    coordinator = await Coordinator.select(
        guild_snowflake=int(source.guild.id),
        member_snowflake=int(member_snowflake),
    )
    if not coordinator:
        raise NotCoordinator
    return True


def coordinator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_sysadmin_wrapper,
            is_developer_wrapper,
            is_guild_owner_wrapper,
            is_administrator_wrapper,
            is_coordinator_at_all,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "Member is not a sysadmin, developer, guild owner, administrator or coordinator in this channel."
        )

    predicate._permission_level = "Coordinator"
    return commands.check(predicate)


class CoordinatorService:

    lines, pages = [], []

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        coordinators = await Coordinator.select(singular=False, **where_kwargs)
        for coordinator in coordinators:
            dictionary.setdefault(coordinator.guild_snowflake, {"members": {}})
            dictionary[coordinator.guild_snowflake]["members"].setdefault(
                coordinator.member_snowflake, {"coordinators": {}}
            )
            dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"].setdefault(coordinator.channel_snowflake, {})
            dictionary[coordinator.guild_snowflake]["members"][
                coordinator.member_snowflake
            ]["coordinators"][coordinator.channel_snowflake].update(
                {"placeholder": "placeholder"}
            )
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                CoordinatorService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                CoordinatorService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Coordinators {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await CoordinatorService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, coordinator_dictionary in guild_data.get(
                "members"
            ).items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    CoordinatorService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in coordinator_dictionary.get(
                    "coordinators"
                ).items():
                    if not isinstance(
                        object_dict.get("object", None), discord.abc.GuildChannel
                    ):
                        channel = guild.get_channel(channel_snowflake)
                        CoordinatorService.lines.append(
                            f"**Channel:** {channel.mention}"
                        )
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(CoordinatorService.lines),
                            inline=False,
                        )
                        embed = flush_page(
                            embed, CoordinatorService.pages, title, guild.name
                        )
                        CoordinatorService.lines = []
                        field_count = 0
            if CoordinatorService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(CoordinatorService.lines),
                    inline=False,
                )
            CoordinatorService.pages.append(embed)
        return CoordinatorService.pages

    @classmethod
    async def toggle_coordinator(cls, channel_dict, member_dict, snowflake_kwargs):
        from vyrtuous.utils.check import has_equal_or_lower_role

        await has_equal_or_lower_role(
            snowflake_kwargs=snowflake_kwargs,
            member_snowflake=member_dict.get("id", None),
        )
        where_kwargs = {}
        where_kwargs.update(channel_dict.get("columns", None))
        where_kwargs.update(member_dict.get("columns", None))
        coordinator = await Coordinator.select(singular=True, **where_kwargs)
        if coordinator:
            await Coordinator.delete(**where_kwargs)
            action = "revoked"
        else:
            coordinator = Coordinator(**where_kwargs)
            await coordinator.create()
            action = "granted"
        return (
            f"Coordinator access has been {action} for {member_dict.get("mention", None)} "
            f"in {channel_dict.get("mention", None)}."
        )
