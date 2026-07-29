"""!/bin/python3
administrator_role_service.py The purpose of this program is to service administrator roles.

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

# import discord
# from discord.ext import commands
#
# from vyrtuous.bot.discord_bot import DiscordBot
# from vyrtuous.db.administrator import AdministratorRole
# from vyrtuous.db.database_factory import DatabaseFactory
# from vyrtuous.listing import list_service
# from vyrtuous.utils.errors.error import GuildNotFound, MemberNotFound, RoleNotFound
# from vyrtuous.utils.messaging import emojis
# from vyrtuous.utils.tracking import data_builder, stream_service
# from vyrtuous.utils.users import administrator_service
#
# MODEL = AdministratorRole
#
#
# async def is_added_role_administrator(
#     guild_snowflake: int, role_snowflake: int
# ) -> bool:
#     database_factory: DatabaseFactory = DatabaseFactory(MODEL)
#     where_kwargs = {
#         "guild_snowflake": int(guild_snowflake),
#         "role_snowflake": int(role_snowflake),
#     }
#     administrator_roles = await database_factory.select(
#         singular=False,
#         **where_kwargs,
#     )
#     if not administrator_roles:
#         return False
#     return True
#
#
# async def toggle_administrator_role(
#     author_snowflake: int,
#     guild_snowflake: int,
#     message_snowflake: int,
#     message_channel_snowflake: int,
#     role_snowflake: int,
# ) -> list[discord.Embed]:
#     bot: DiscordBot = DiscordBot.get_instance()
#     guild = bot.get_guild(guild_snowflake)
#     if guild is None:
#         raise GuildNotFound(str(guild_snowflake))
#     role = guild.get_role(role_snowflake)
#     if role is None:
#         raise RoleNotFound(str(role_snowflake))
#     database_factory: DatabaseFactory = DatabaseFactory(MODEL)
#     title = f"{emojis.get_random_emoji()} Administrators and Roles"
#     administrators = await administrator_service.administrators_by_role(
#         role_snowflake=role_snowflake
#     )
#     administrator_roles = await is_added_role_administrator(
#         guild_snowflake=guild_snowflake, role_snowflake=role_snowflake
#     )
#     if administrator_roles:
#         action = "revoked"
#         if administrator_roles:
#             await database_factory.delete(role_snowflake=role_snowflake)
#         revoked_members: dict[int, dict[int, list[discord.Member]]] = {}
#         for administrator in administrators:
#             member = role.guild.get_member(administrator.member_snowflake)
#             if member is None:
#                 raise MemberNotFound(str(administrator.member_snowflake))
#             await administrator_service.removed_role(
#                 guild_snowflake=guild_snowflake,
#                 member_snowflake=administrator.member_snowflake,
#                 role_snowflake=role_snowflake,
#             )
#             await log_xadmin(
#                 author_snowflake=author_snowflake,
#                 display=True,
#                 guild_snowflake=role.guild.id,
#                 member_snowflake=member.id,
#                 message_snowflake=message_snowflake,
#                 message_channel_snowflake=message_channel_snowflake,
#                 role_snowflake=role_snowflake,
#             )
#
#             revoked_members.setdefault(guild_snowflake, {}).setdefault(
#                 role_snowflake, []
#             ).append(member)
#         members = revoked_members.get(guild_snowflake, {}).get(role_snowflake, [])
#     else:
#         action = "granted"
#         granted_members: dict[int, dict[int, list[discord.Member]]] = {}
#         granted_members.setdefault(role.guild.id, {})[role_snowflake] = []
#         administrator_role = AdministratorRole(
#             guild_snowflake=role.guild.id, role_snowflake=role_snowflake
#         )
#         await database_factory.create(administrator_role)
#         for member in role.members:
#             await administrator_service.added_role(
#                 guild_snowflake=role.guild.id,
#                 member_snowflake=member.id,
#                 role_snowflake=role_snowflake,
#             )
#             await log_admin(
#                 author_snowflake=author_snowflake,
#                 display=True,
#                 guild_snowflake=role.guild.id,
#                 member_snowflake=member.id,
#                 message_snowflake=message_snowflake,
#                 message_channel_snowflake=message_channel_snowflake,
#                 role_snowflake=role_snowflake,
#             )
#             granted_members[role.guild.id][role_snowflake].append(member)
#         members = granted_members.get(role.guild.id, {}).get(role_snowflake, [])
#     embed = discord.Embed(
#         title=title,
#         description=f"`{role.name}` was `{action}`.",
#         color=discord.Color.red() if action == "revoked" else discord.Color.green(),
#     )
#     embed.add_field(name="Role ID", value=str(role_snowflake), inline=False)
#     embed.add_field(name="Guild", value=str(guild.name), inline=False)
#
#     chunks = []
#     chunk = []
#     pages: list[discord.Embed] = []
#     for member in members:
#         chunk.append(member)
#         if len(chunk) == list_service.CHUNK_SIZE:
#             chunks.append(chunk)
#             chunk = []
#     if chunk:
#         chunks.append(chunk)
#     field_count = 1
#     page_number = 1
#     for chunk in chunks:
#         embed = discord.Embed(
#             title=f"Members {action.capitalize()}",
#             color=(
#                 discord.Color.red() if action == "revoked" else discord.Color.green()
#             ),
#         )
#         for member in chunk:
#             embed.add_field(
#                 name=f"{field_count}. {member}",
#                 value=f"{member.mention} ({member.id})",
#                 inline=False,
#             )
#             field_count += 1
#         embed.set_footer(text=f"Page {page_number}")
#         pages.append(embed)
#         page_number += 1
#     return pages
#
#
# async def get_administrator_roles(guild_snowflake=None) -> list[int]:
#     database_factory: DatabaseFactory = DatabaseFactory(MODEL)
#     roles = await database_factory.select(
#         guild_snowflake=guild_snowflake, singular=False
#     )
#     role_snowflakes = []
#     for role in roles:
#         role_snowflakes.append(role.role_snowflake)
#     return role_snowflakes
#
#
# async def log_admin(
#     author_snowflake: int | None,
#     display: bool,
#     guild_snowflake: int,
#     member_snowflake: int,
#     message_snowflake: int | None,
#     message_channel_snowflake: int | None,
#     role_snowflake: int,
# ):
#     channel_snowflake = None
#     duration = None
#     is_channel_scope = None
#     reason = None
#     target = None
#     await data_builder.save_data(
#         author_snowflake=author_snowflake or None,
#         channel_snowflake=channel_snowflake,
#         duration=duration or None,
#         guild_snowflake=guild_snowflake,
#         identifier="admin",
#         member_snowflake=member_snowflake,
#         reason=reason or "No reason provided.",
#         role_snowflake=role_snowflake or None,
#         target=target or None,
#     )
#     if display:
#         await stream_service.send_log(
#             author_snowflake=author_snowflake or None,
#             channel_snowflake=channel_snowflake,
#             identifier="admin",
#             duration=duration or None,
#             guild_snowflake=guild_snowflake,
#             is_channel_scope=is_channel_scope,
#             member_snowflake=member_snowflake,
#             message_snowflake=message_snowflake or None,
#             message_channel_snowflake=message_channel_snowflake or None,
#             reason=reason or "No reason provided.",
#             role_snowflake=role_snowflake or None,
#             target=target or None,
#         )
#
#
# async def log_xadmin(
#     author_snowflake: int | None,
#     display: bool,
#     guild_snowflake: int,
#     member_snowflake: int,
#     message_snowflake: int | None,
#     message_channel_snowflake: int | None,
#     role_snowflake: int,
# ):
#     channel_snowflake = None
#     duration = None
#     is_channel_scope = None
#     reason = None
#     target = None
#     await data_builder.save_data(
#         author_snowflake=author_snowflake or None,
#         channel_snowflake=channel_snowflake,
#         duration=duration or None,
#         guild_snowflake=guild_snowflake,
#         identifier="xadmin",
#         member_snowflake=member_snowflake,
#         reason=reason or "No reason provided.",
#         role_snowflake=role_snowflake or None,
#         target=target or None,
#     )
#     if display:
#         await stream_service.send_log(
#             author_snowflake=author_snowflake or None,
#             channel_snowflake=channel_snowflake,
#             identifier="xadmin",
#             duration=duration or None,
#             guild_snowflake=guild_snowflake,
#             is_channel_scope=is_channel_scope,
#             member_snowflake=member_snowflake,
#             message_snowflake=message_snowflake or None,
#             message_channel_snowflake=message_channel_snowflake or None,
#             reason=reason or "No reason provided.",
#             role_snowflake=role_snowflake or None,
#             target=target or None,
#         )
