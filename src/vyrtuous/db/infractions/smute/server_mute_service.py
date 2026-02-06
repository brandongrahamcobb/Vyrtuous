"""server_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the server mute moderation.

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

from vyrtuous.base.service import Service
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.infractions.smute.server_mute import ServerMute
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


class ServerMuteService(Service):

    lines, pages = [], []

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        server_mutes = await ServerMute.select(singular=False, **where_kwargs)
        for server_mute in server_mutes:
            dictionary.setdefault(server_mute.guild_snowflake, {"members": {}})
            dictionary[server_mute.guild_snowflake]["members"].setdefault(
                server_mute.member_snowflake, {"server_mutes": {}}
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
                ServerMuteService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                ServerMuteService.pages.extend(
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
        title = f"{get_random_emoji()} Server Mutes {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get("object", None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await ServerMuteService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        smute_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, server_mute_dictionary in guild_data.get(
                "members", {}
            ).items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    ServerMuteService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                    field_count += 1
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                smute_n += 1
                field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(ServerMuteService.lines),
                        inline=False,
                    )
                    embed = flush_page(
                        embed, ServerMuteService.pages, title, guild.name
                    )
                    ServerMuteService.lines = []
                    field_count = 0
            if ServerMuteService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(ServerMuteService.lines),
                    inline=False,
                )
            ServerMuteService.pages.append(embed)
        if ServerMuteService.pages:
            ServerMuteService.pages[0].description = f'**({smute_n})**'
        return ServerMuteService.pages

    @classmethod
    async def toggle_server_mute(cls, default_kwargs, member_dict, reason):
        bot = DiscordBot.get_instance()
        updated_kwargs = default_kwargs.copy()
        guild_snowflake = default_kwargs.get("guild_snowflake", None)
        guild = bot.get_guild(guild_snowflake)
        await PermissionService.has_equal_or_lower_role(
            updated_kwargs=default_kwargs,
            member_snowflake=member_dict.get("id", None),
        )
        updated_kwargs.update(member_dict.get("columns", None))
        del updated_kwargs["channel_snowflake"]
        server_mute = await ServerMute.select(singular=True, **updated_kwargs)
        if not server_mute:
            server_mute = ServerMute(
                **updated_kwargs,
                reason=reason,
            )
            await server_mute.create()
            action = "muted"
            should_be_muted = True
        else:
            await ServerMute.delete(**updated_kwargs)
            action = "unmuted"
            should_be_muted = False

        if (
            member_dict.get("object", None).voice
            and member_dict.get("object", None).voice.channel
        ):
            await member_dict.get("object", None).edit(mute=should_be_muted)
        return f"Successfully server {action} {member_dict.get("mention", None)} in {guild.name}."
