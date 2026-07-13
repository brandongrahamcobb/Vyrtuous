"""!/bin/python3
developer_service.py The purpose of this program is to extend Service to service the developer class.

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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.developer import Developer, NotDeveloper
from vyrtuous.utils.tracking import data_builder

# from vyrtuous.models.duration import DurationBuilder
# from vyrtuous.utils.tracking import bug_service

MODEL = Developer


async def is_developer(member_snowflake: int) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    developer = await database_factory.select(
        member_snowflake=int(member_snowflake), singular=True
    )
    if not developer:
        raise NotDeveloper
    return True


async def report_issue(author, message, reference) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    online_developer_mentions = []
    sysadmin = bot.get_user(bot.config["discord_owner_id"])
    if sysadmin:
        online_developer_mentions.append(sysadmin.mention)
    if message.guild:
        msg = f"Issue reported by {author.name}!\n**Message:** {message.jump_url}\n**Reference:** {reference}"
        developers = await database_factory.select(singular=False)
        for developer in developers:
            member = message.guild.get_member(developer.member_snowflake)
            if member and member.status != discord.Status.offline:
                online_developer_mentions.append(member.mention)
                try:
                    await member.send(msg)
                    if sysadmin:
                        await sysadmin.send(msg)
                except discord.Forbidden as e:
                    bot.logger.warning(
                        f"Unable to send a developer log ID: {id}. {str(e).capitalize()}"
                    )
    msg = "Your report has been submitted"
    if online_developer_mentions:
        msg = f"{message}. The developers {', '.join(online_developer_mentions)} are online and will respond to your report shortly."
    await author.send(msg)


async def toggle_developer(
    author_snowflake: int,
    member_snowflake: int,
    message_snowflake: int,
    message_channel_snowflake: int,
) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    found = False
    developers = await database_factory.select(singular=False)
    found = any(
        developer.member_snowflake == member_snowflake for developer in developers
    )
    member = bot.get_user(member_snowflake)
    if member:
        member_str = member.mention
    else:
        simplified_member = bot.registry.get(MemberState).active.get(
            member_snowflake, None
        )
        if simplified_member:
            member_str = simplified_member[0]
        else:
            raise commands.MemberNotFound(str(member_snowflake))
    if found:
        await database_factory.delete(member_snowflake=member_snowflake)
        await log_xdev(
            author_snowflake=author_snowflake,
            member_snowflake=member_snowflake,
        )
        action = "revoked"
    else:
        new_developer = MODEL(member_snowflake=member_snowflake)
        await database_factory.create(new_developer)
        await log_dev(
            author_snowflake=author_snowflake,
            member_snowflake=member_snowflake,
        )
        action = "granted"
    return f"Developer access for {member_str} has been {action} globally."


async def log_dev(
    author_snowflake: int | None,
    member_snowflake: int,
):
    channel_snowflake = None
    duration_value = None
    guild_snowflake = None
    reason = None
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake or None,
        channel_snowflake=channel_snowflake,
        duration_value=duration_value or None,
        guild_snowflake=guild_snowflake,
        identifier="dev",
        member_snowflake=member_snowflake,
        reason=reason or "No reason provided.",
        role_snowflake=role_snowflake or None,
        target=target or None,
    )


async def log_xdev(
    author_snowflake: int | None,
    member_snowflake: int,
):
    channel_snowflake = None
    duration_value = None
    guild_snowflake = None
    reason = None
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake or None,
        channel_snowflake=channel_snowflake,
        duration_value=duration_value or None,
        guild_snowflake=guild_snowflake,
        identifier="xdev",
        member_snowflake=member_snowflake,
        reason=reason or "No reason provided.",
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
