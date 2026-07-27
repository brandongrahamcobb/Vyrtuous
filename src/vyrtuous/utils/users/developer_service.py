"""!/bin/python3
developer_service.py The purpose of this program is to extend Service to service the developer class.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.developer import Developer, NotDeveloper
from vyrtuous.models.duration import DurationObject
from vyrtuous.utils.tracking import data_builder, stream_service

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
            display=True,
            member_snowflake=member_snowflake,
            message_channel_snowflake=message_channel_snowflake,
            message_snowflake=message_snowflake,
        )
        action = "revoked"
    else:
        new_developer = MODEL(member_snowflake=member_snowflake)
        await database_factory.create(new_developer)
        await log_dev(
            author_snowflake=author_snowflake,
            display=True,
            member_snowflake=member_snowflake,
            message_channel_snowflake=message_channel_snowflake,
            message_snowflake=message_snowflake,
        )
        action = "granted"
    return f"Developer access for {member_str} has been {action} globally."


async def log_dev(
    author_snowflake: int,
    display: bool,
    member_snowflake: int,
    message_channel_snowflake: int,
    message_snowflake: int,
):
    channel_snowflake = None
    duration = DurationObject(number=0, prefix="", sign=1, unit="")
    guild_snowflake = None
    is_channel_scope = False
    reason = "No reason provided."
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake,
        channel_snowflake=channel_snowflake or None,
        duration=duration,
        guild_snowflake=guild_snowflake or None,
        identifier="dev",
        member_snowflake=member_snowflake,
        reason=reason,
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake,
            channel_snowflake=channel_snowflake,
            identifier="dev",
            duration=duration,
            guild_snowflake=guild_snowflake,
            is_channel_scope=is_channel_scope,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake or None,
            message_channel_snowflake=message_channel_snowflake or None,
            reason=reason,
            role_snowflake=role_snowflake or None,
            target=target or None,
        )


async def log_xdev(
    author_snowflake: int,
    display: bool,
    member_snowflake: int,
    message_channel_snowflake: int,
    message_snowflake: int,
):
    channel_snowflake = None
    duration = DurationObject(number=0, prefix="", sign=1, unit="")
    guild_snowflake = None
    is_channel_scope = False
    reason = "No reason provided."
    role_snowflake = None
    target = None
    await data_builder.save_data(
        author_snowflake=author_snowflake,
        channel_snowflake=channel_snowflake or None,
        duration=duration,
        guild_snowflake=guild_snowflake or None,
        identifier="xdev",
        member_snowflake=member_snowflake,
        reason=reason,
        role_snowflake=role_snowflake or None,
        target=target or None,
    )
    if display:
        await stream_service.send_log(
            author_snowflake=author_snowflake,
            channel_snowflake=channel_snowflake,
            identifier="xdev",
            duration=duration,
            guild_snowflake=guild_snowflake,
            is_channel_scope=is_channel_scope,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake or None,
            message_channel_snowflake=message_channel_snowflake or None,
            reason=reason,
            role_snowflake=role_snowflake or None,
            target=target or None,
        )
