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
from vyrtuous.models.duration import DurationBuilder
from vyrtuous.utils.tracking import bug_service

MODEL = Developer


async def is_developer(member_snowflake: int) -> bool:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    developer = await database_factory.select(
        member_snowflake=int(member_snowflake), singular=True
    )
    if not developer:
        raise NotDeveloper
    return True


async def is_developer_wrapper(context) -> bool:
    return await is_developer(member_snowflake=int(context.member_snowflake))


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


async def toggle_developer(member_snowflake: int) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    found = False
    developers = await database_factory.select(singular=False)
    found = any(
        developer.member_snowflake == member_snowflake for developer in developers
    )
    if found:
        await database_factory.delete(member_snowflake=member_snowflake)
        action = "revoked"
    else:
        new_developer = MODEL(member_snowflake=member_snowflake)
        await database_factory.create(new_developer)
        action = "granted"
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
    return f"Developer access for {member_str} has been {action} globally."


async def ping_about_expired_bugs(
    channel,
    embed,
    member,
    member_snowflakes,
    msg,
    notes,
    updated_at,
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    duration_builder = DurationBuilder()
    assigned_developer_mentions = []
    for developer_snowflake in member_snowflakes:
        assigned_developer = bot.get_user(developer_snowflake)
        if assigned_developer:
            assigned_developer_mentions.append(assigned_developer.mention)
    embed.add_field(
        name=f"Updated: {duration_builder.from_timestamp(updated_at).to_unix_ts()}",
        value=f"**Link:** {msg.jump_url}\n**Developers:** {', '.join(assigned_developer_mentions)}\n**Notes:** {notes}",
        inline=False,
    )
    for developer_snowflake in member_snowflakes:
        member = bot.get_user(developer_snowflake)
        if member is None:
            bot.logger.info(
                f"Unable to locate member {developer_snowflake} in channel {channel.name} ({channel.id}) not sending developer log."
            )
            continue
        try:
            await member.send(embed=embed)
            bot.logger.info(
                f"Sent the issue to member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id})."
            )
        except Exception as e:
            bot.logger.warning(
                f"Unable to send the issue to member {member.display_name} ({member.id}) in channel {channel.name} ({channel.id}). {str(e).capitalize()}"
            )


async def handle_developer_assignment(member, reference) -> discord.Embed:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    developers = await database_factory.select(singular=False)
    for developer in developers:
        if developer.member_snowflake == member.id:
            bug, state = await bug_service.handle_bug_assignment(
                developer=developer, reference=reference
            )
            if not state:
                action = "unassigned"
            else:
                action = "assigned"
            return await bug_service.create_embed(
                action=action,
                bug=bug,
                member=member,
            )
