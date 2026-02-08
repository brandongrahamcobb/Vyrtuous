"""!/bin/python3
role_service.py The purpose of this program is to extend AliasService to service the role class.

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

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.base.record_service import RecordService
from vyrtuous.db.mgmt.stream.stream_service import StreamService
from vyrtuous.db.roles.role import Role
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class RoleService(RecordService):
    lines, pages = [], []
    model = Role

    @classmethod
    async def act_embed(cls, ctx):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(ctx.target_channel_snowflake)
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        role = guild.get_role(ctx.target_role_snowflake)
        embed = discord.Embed(
            title=f"{get_random_emoji()} {member.display_name} has been granted a role",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, ctx):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(ctx.target_channel_snowflake)
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        role = guild.get_role(ctx.target_role_snowflake)
        embed = discord.Embed(
            title=f"{get_random_emoji()} {member.display_name}'s role has been revoked",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def administer_role(cls, guild_snowflake, member_snowflake, role_snowflake):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        role = guild.get_role(role_snowflake)
        try:
            await member.add_roles(role, reason="Granting role.")
        except discord.Forbidden as e:
            logger.error(str(e).capitalize())

    @classmethod
    async def revoke_role(cls, guild_snowflake, member_snowflake, role_snowflake):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        role = guild.get_role(role_snowflake)
        try:
            await member.remove_roles(role, reason="Revoking role.")
        except discord.Forbidden as e:
            logger.error(str(e).capitalize())

    @classmethod
    async def added_role(
        cls,
        category_class,
        category_role_class,
        guild_snowflake,
        member_snowflake,
        role_snowflake,
    ):
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": str(role_snowflake),
        }
        role = await category_role_class.select(singular=True, **kwargs)
        if role:
            kwargs.update({"channel_snowflake": role.channel_snowflake})
            kwargs.update({"member_snowflake": role.member_snowflake})
            msg = f"Member ({member_snowflake}) was granted the role ({role_snowflake}) for category ({category_class.__name__()}) related to channel ({role.channel_snowflake}) in guild ({guild_snowflake})."
            category = category_class(**kwargs)
            await category.create()
            logger.info(msg)
        else:
            return

    @classmethod
    async def removed_role(
        cls,
        category_class,
        category_role_class,
        guild_snowflake,
        member_snowflake,
        role_snowflake,
    ):
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "role_snowflake": str(role_snowflake),
        }
        role = await category_role_class.select(singular=True, **kwargs)
        if role:
            kwargs.update({"channel_snowflake": role.channel_snowflake})
            kwargs.update({"member_snowflake": role.member_snowflake})
            msg = f"Member ({member_snowflake}) was revoked the role ({role_snowflake}) for category ({category_class.__name__()}) related to channel ({role.channel_snowflake}) in guild ({guild_snowflake})."
            await category_class.delete(**kwargs)
            logger.info(msg)
        else:
            return

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        roles = await Role.select(singular=False, **where_kwargs)
        for role in roles:
            dictionary.setdefault(role.guild_snowflake, {"members": {}})
            dictionary[role.guild_snowflake]["members"].setdefault(
                role.member_snowflake, {"roles": {}}
            )
            dictionary[role.guild_snowflake]["members"][role.member_snowflake][
                "roles"
            ].setdefault(role.channel_snowflake, {})
            dictionary[role.guild_snowflake]["members"][role.member_snowflake]["roles"][
                role.channel_snowflake
            ] = {"role": role.role_snowflake}
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                RoleService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                RoleService.pages.extend(
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
        title = f"{get_random_emoji()} Role {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await Role.build_cleaned_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        role_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, role_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    RoleService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in role_dictionary.get(
                    "roles"
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    role = channel_dictionary.get("role", None)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        RoleService.lines.append(f"**Channel:** {channel.mention}")
                    RoleService.lines.append(f"**Role:** {role.mention}")
                    role_n += 1
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(RoleService.lines),
                            inline=False,
                        )
                        embed = flush_page(embed, RoleService.pages, title, guild.name)
                        RoleService.lines = []
                        field_count = 0
            if RoleService.lines:
                embed.add_field(
                    name="Information", value="\n".join(RoleService.lines), inline=False
                )
            RoleService.pages.append(embed)
        if RoleService.pages:
            RoleService.pages[0].description = f"**({role_n})**"
        return RoleService.pages

    @classmethod
    async def enforce(cls, ctx, source, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        added_role = Role(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            role_snowflake=ctx.target_role_snowflake,
        )
        await added_role.create()
        role = guild.get_role(ctx.target_role_snowflake)
        if role:
            await RoleService.administer_role(
                guild_snowflake=ctx.source_guild_snowflake,
                member_snowflake=ctx.target_member_snowflake,
                role_snowflake=ctx.target_role_snowflake,
            )
        await StreamService.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="role",
            member=member,
            source=source,
        )
        embed = await RoleService.act_embed(ctx=ctx)
        return await state.end(success=embed)

    @classmethod
    async def undo(cls, ctx, source, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        await Role.delete(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            role_snowflake=ctx.target_role_snowflake,
        )
        role = guild.get_role(ctx.target_role_snowflake)
        if role:
            await RoleService.revoke_role(
                guild_snowflake=ctx.source_guild_snowflake,
                member_snowflake=ctx.target_member_snowflake,
                role_snowflake=ctx.target_role_snowflake,
            )
        await StreamService.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="unrole",
            is_modification=True,
            member=member,
            source=source,
        )
        embed = await RoleService.undo_embed(ctx=ctx)
        return await state.end(success=embed)
