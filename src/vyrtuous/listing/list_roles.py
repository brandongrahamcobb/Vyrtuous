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

# from dataclasses import dataclass, field
# from typing import Dict, List
#
# import discord
#
# from vyrtuous.db.database_factory import DatabaseFactory
# from vyrtuous.bot.discord_bot import DiscordBot
# from vyrtuous.listing import list_service
#
# MODEL = Role
#
#
# @dataclass
# class RoleDictionary:
#     data: Dict[int, Dict[str, Dict[int, Dict[str, Dict[int, Dict[str, int]]]]]] = field(
#         default_factory=dict
#     )
#     skipped_guilds: List[discord.Embed] = field(default_factory=list)
#     skipped_members: List[discord.Embed] = field(default_factory=list)
#
#
# async def build_dictionary(obj):
#     database_factory = DatabaseFactory(MODEL)
#     roles = []
#     dictionary = {}
#     if isinstance(obj, discord.Guild):
#         roles = await database_factory.select(guild_snowflake=obj.id, singular=False)
#     elif isinstance(obj, discord.abc.GuildChannel):
#         roles = await database_factory.select(channel_snowflake=obj.id, singular=False)
#     elif isinstance(obj, discord.Member):
#         roles = await database_factory.select(member_snowflake=obj.id, singular=False)
#     else:
#         roles = await database_factory.select(singular=False)
#     if roles:
#         for role in roles:
#             dictionary.setdefault(role.guild_snowflake, {"members": {}})
#             dictionary[role.guild_snowflake]["members"].setdefault(
#                 role.member_snowflake, {"roles": {}}
#             )
#             dictionary[role.guild_snowflake]["members"][role.member_snowflake][
#                 "roles"
#             ].setdefault(role.channel_snowflake, {})
#             dictionary[role.guild_snowflake]["members"][role.member_snowflake]["roles"][
#                 role.channel_snowflake
#             ] = {"role": role.role_snowflake}
#     return dictionary
#
#
# async def build_pages(is_at_home: bool, obj):
#     bot: DiscordBot = DiscordBot.get_instance()
#     lines, pages: list[discord.Embed] = []
#     title = f"{emojis.get_random_emoji()} Role {f'for {obj.name}' if isinstance(obj, discord.Member) else ''}"
#
#     dictionary = await build_dictionary(obj=obj)
#     processed_dictionary = await list_service.process_dictionary(
#         cls=RoleDictionary, dictionary=dictionary
#     )
#
#     for guild_snowflake, guild_data in processed_dictionary.data.items():
#         role_n = 0
#         field_count = 0
#         lines = []
#         thumbnail = False
#         guild = bot.get_guild(guild_snowflake)
#         if guild is None:
#             continue
#         embed = discord.Embed(
#             title=title, description=guild.name, color=discord.Color.blue()
#         )
#         for member_snowflake, role_dictionary in guild_data.get("members").items():
#             member = guild.get_member(member_snowflake)
#             if member is None:
#                 continue
#             if not isinstance(obj, discord.Member):
#                 lines.append(f"**User:** {member.display_name} {member.mention}")
#             elif not thumbnail:
#                 embed.set_thumbnail(url=obj.display_avatar.url)
#                 thumbnail = True
#             for channel_snowflake, channel_dictionary in role_dictionary.get(
#                 "roles"
#             ).items():
#                 channel = guild.get_channel(channel_snowflake)
#                 if channel is None:
#                     continue
#                 role = channel_dictionary.get("role", None)
#                 if not isinstance(obj, discord.abc.GuildChannel):
#                     lines.append(f"**Channel:** {channel.mention}")
#                 lines.append(f"**Role:** {role.mention}")
#                 role_n += 1
#                 field_count += 1
#                 if field_count >= list_service.CHUNK_SIZE:
#                     embed.add_field(
#                         name="Information",
#                         value="\n".join(lines),
#                         inline=False,
#                     )
#                     embed = list_service.flush_page(
#                         embed, pages, title, guild.name
#                     )
#                     lines = []
#                     field_count = 0
#         if lines:
#             embed.add_field(name="Information", value="\n".join(lines), inline=False)
#         original_description = embed.description or ""
#         embed.description = f"**{original_description} ({role_n})**"
#         pages.append(embed)
#     if is_at_home:
#         pages.extend(processed_dictionary.skipped_guilds)
#         pages.extend(processed_dictionary.skipped_members)
#     return pages
