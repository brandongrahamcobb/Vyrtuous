"""moderator.py The purpose of this program is to inherit from the DatabaseFactory to provide the moderator role.

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

from vyrtuous.db.roles.admin.administrator_service import is_administrator_wrapper
from vyrtuous.db.roles.coord.coordinator_service import is_coordinator_at_all
from vyrtuous.db.roles.dev.developer_service import is_developer_wrapper
from vyrtuous.db.roles.mod.moderator import Moderator
from vyrtuous.db.roles.owner.guild_owner_service import is_guild_owner_wrapper
from vyrtuous.db.roles.sysadmin.sysadmin_service import is_sysadmin_wrapper
from vyrtuous.discord.author import resolve_author
from vyrtuous.discord.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import CHUNK_SIZE
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
class NotModerator(commands.CheckFailure):
    def __init__(
        self,
        message="Member is not a sysadmin, developer, guild owner, administrator, coordinator or moderator in this channel.",
    ):
        super().__init__(message)


async def is_moderator_wrapper(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    member = resolve_author(source=source)
    member_snowflake = member.id
    return await is_moderator(
        channel_snowflake=source.channel.id,
        guild_snowflake=source.guild.id,
        member_snowflake=int(member_snowflake),
    )


async def is_moderator(
    channel_snowflake: int, guild_snowflake: int, member_snowflake: int
) -> bool:
    moderator = await Moderator.select(
        channel_snowflake=int(channel_snowflake),
        guild_snowflake=int(guild_snowflake),
        member_snowflake=int(member_snowflake),
        singular=True,
    )
    if not moderator:
        raise NotModerator
    return True


async def is_moderator_at_all(
    source: Union[commands.Context, discord.Interaction, discord.Message],
) -> bool:
    member = resolve_author(source=source)
    member_snowflake = member.id
    moderator = await Moderator.select(
        guild_snowflake=int(source.guild.id), member_snowflake=int(member_snowflake)
    )
    if not moderator:
        raise NotModerator
    return True


def moderator_predicator():
    async def predicate(
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        for verify in (
            is_sysadmin_wrapper,
            is_developer_wrapper,
            is_guild_owner_wrapper,
            is_administrator_wrapper,
            is_coordinator_at_all,
            is_moderator_at_all,
        ):
            try:
                if await verify(source):
                    return True
            except commands.CheckFailure:
                continue
        raise commands.CheckFailure(
            "Member is not a sysadmin, developer, guild owner, administrator, coordinator or moderator in this channel."
        )

    predicate._permission_level = "Moderator"
    return commands.check(predicate)


class ModeratorService:

    lines, pages = [], []

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        moderators = await Moderator.select(singular=False, **where_kwargs)
        for moderator in moderators:
            dictionary.setdefault(moderator.guild_snowflake, {"members": {}})
            dictionary[moderator.guild_snowflake]["members"].setdefault(
                moderator.member_snowflake, {"moderators": {}}
            )
            dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"].setdefault(moderator.channel_snowflake, {})
            dictionary[moderator.guild_snowflake]["members"][
                moderator.member_snowflake
            ]["moderators"][moderator.channel_snowflake].update(
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
                ModeratorService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                ModeratorService.pages.extend(
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
        title = f"{get_random_emoji()} Moderator {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await ModeratorService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, member_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not isinstance(object_dict.get("object", None), discord.Member):
                    ModeratorService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in member_dictionary.get(
                    "moderators", {}
                ).items():
                    if not isinstance(
                        object_dict.get("object", None), discord.abc.GuildChannel
                    ):
                        channel = guild.get_channel(channel_snowflake)
                        ModeratorService.lines.append(f"**Channel:** {channel.mention}")
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(ModeratorService.lines),
                            inline=False,
                        )
                        embed = flush_page(
                            embed, ModeratorService.pages, title, guild.name
                        )
                        ModeratorService.lines = []
                        field_count = 0
            if ModeratorService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(ModeratorService.lines),
                    inline=False,
                )
            ModeratorService.pages.append(embed)
        return ModeratorService.pages

    @classmethod
    async def toggle_moderator(cls, channel_dict, member_dict, snowflake_kwargs):
        from vyrtuous.permissions.check import has_equal_or_lower_role

        await has_equal_or_lower_role(
            snowflake_kwargs=snowflake_kwargs,
            member_snowflake=member_dict.get("id", None),
        )
        where_kwargs = {}
        where_kwargs.update(channel_dict.get("columns", None))
        where_kwargs.update(member_dict.get("columns", None))
        moderator = await Moderator.select(singular=True, **where_kwargs)
        if moderator:
            await Moderator.delete(**where_kwargs)
            action = "revoked"
        else:
            moderator = Moderator(**where_kwargs)
            await moderator.create()
            action = "granted"
        return (
            f"Moderator access for {member_dict.get("mention", None)} has been "
            f"{action} in {channel_dict.get("mention", None)}."
        )
