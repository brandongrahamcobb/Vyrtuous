"""!/bin/python3
bug_service.py The purpose of this program is to extend Service to service the bug command class.

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
from vyrtuous.db.bug import Bug
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.messaging import emojis

MODEL = Bug


async def update_bug(action: str, notes, reference) -> str:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    message = "You successfully "
    bug = await database_factory.select(id=reference, resolved=False, singular=True)
    if not bug:
        return f"Unresolved issue not found for reference ({reference})."
    if action and action.lower() == "resolve":
        where_kwargs = {"id": reference}
        set_kwargs = {"resolved": True}
        await database_factory.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
        message += "resolved the issue. The record will remain in the database for the next 30 days."
    elif action and action.lower() == "append":
        where_kwargs = {"id": reference}
        set_kwargs = {"notes": bug.notes + notes if bug.notes else notes}
        await database_factory.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
        message += "appended to the previous notes."
    elif action and action.lower() == "overwrite":
        where_kwargs = {"id": reference}
        set_kwargs = {"notes": notes}
        await database_factory.update(where_kwargs=where_kwargs, set_kwargs=set_kwargs)
        message += "overwrote the previous notes."
    return message


async def create_embed(action: str, bug, member) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    current_developer_mentions = []
    guild = bot.get_guild(bug.guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(bug.guild_snowflake)
    for current_developer_snowflake in bug.member_snowflakes:
        current_developer = guild.get_member(current_developer_snowflake)
        if current_developer is None:
            continue
        current_developer_mentions.append(current_developer.mention)
    channel = guild.get_channel(bug.channel_snowflake)
    if channel is None:
        raise commands.ChannelNotFound(bug.channel_snowflake)
    if not isinstance(channel, discord.TextChannel) or not isinstance(
        channel, discord.VoiceChannel
    ):
        return
    toggled_developer = guild.get_member(member.id)
    if toggled_developer is None:
        raise commands.MemberNotFound(member.id)
    try:
        msg = await channel.fetch_message(bug.message_snowflake)
    except discord.NotFound:
        bot.logger.warning(f"Message reference not found ({bug.message_snowflake}).")
    else:
        embed = discord.Embed(
            title=f"{emojis.get_random_emoji()} {toggled_developer.display_name} has been {action}",
            description=(
                f"**Guild:** {guild.name}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Message:** {msg.jump_url}\n"
                f"**Assigned devs:** {', '.join(current_developer_mentions)}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=toggled_developer.display_avatar.url)
        return embed


async def create_bug(message, reference) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    try:
        bug = MODEL(
            channel_snowflake=message.channel.id,
            member_snowflakes=[],
            guild_snowflake=message.guild.id,
            id=reference,
            message_snowflake=message.id,
        )
        await database_factory.create(bug)
    except discord.Forbidden as e:
        bot.logger.info(str(e).capitalize())


async def handle_bug_assignment(developer, reference) -> tuple[Bug, bool]:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    bugs = await database_factory.select(singular=False)
    state = False
    bug = next((bug for bug in bugs if bug.id == reference and not bug.resolved), None)
    if bug:
        where_kwargs = {"id": bug.id}
        member_snowflakes = bug.member_snowflakes
        if developer.member_snowflake in bug.member_snowflakes:
            member_snowflakes.remove(developer.member_snowflake)
            set_kwargs = {"member_snowflakes": member_snowflakes}
            await database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )
            state = False
        else:
            member_snowflakes.append(developer.member_snowflake)
            set_kwargs = {"member_snowflakes": member_snowflakes}
            await database_factory.update(
                set_kwargs=set_kwargs, where_kwargs=where_kwargs
            )
    return bug, state
