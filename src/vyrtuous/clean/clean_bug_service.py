"""!/bin/python3
bug_service.py The purpose of this program is to extend Service to service the bug command class.

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

from datetime import datetime, timedelta, timezone

import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.bug import Bug
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.users import developer_service

MODEL = Bug


async def clean_expired_bugs() -> int:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    count: int = 0
    now = datetime.now(timezone.utc)
    bugs = await database_factory.select(resolved=True, singular=False)
    if bugs:
        for bug in bugs:
            channel_snowflake = int(bug.channel_snowflake)
            guild_snowflake = int(bug.guild_snowflake)
            member_snowflakes = bug.member_snowflakes
            message_snowflake = int(bug.message_snowflake)
            reference = bug.id
            if bug.created_at < now - timedelta(weeks=1):
                guild = bot.get_guild(bug.guild_snowflake)
                if guild is None:
                    bot.logger.info(
                        f"Unable to locate guild {guild_snowflake}, not sending developer log."
                    )
                    continue
                embed = discord.Embed(
                    title=f"\U000026a0\U0000fe0f An issue is unresolved in {guild.name}",
                    color=discord.Color.red(),
                )
                channel = guild.get_channel(channel_snowflake)
                if (
                    channel is None
                    or not isinstance(channel, discord.TextChannel)
                    or not isinstance(channel, discord.VoiceChannel)
                ):
                    bot.logger.info(
                        f"Unable to locate channel {channel_snowflake} in guild {guild.name} ({guild_snowflake}, not sending developer log."
                    )
                    continue
                for member_snowflake in member_snowflakes:
                    member = guild.get_member(member_snowflake)
                    if member is None:
                        bot.logger.info(
                            f"Unable to locate member {member_snowflake} in channel {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), not sending developer log."
                        )
                        continue
                    embed.set_thumbnail(url=member.display_avatar.url)
                    try:
                        msg = await channel.fetch_message(message_snowflake)
                        await developer_service.ping_about_expired_bugs(
                            channel=channel,
                            embed=embed,
                            member=member,
                            member_snowflakes=bug.member_snowflakes,
                            msg=msg,
                            notes=bug.notes,
                            updated_at=bug.updated_at,
                        )
                    except Exception as e:
                        bot.logger.warning(
                            f"Unable to locate a message {message_snowflake} in {channel.name} ({channel.id}) in guild {guild.name} ({guild_snowflake}), deleting developer log. {str(e).capitalize()}"
                        )
                    await database_factory.delete(id=reference)
                    count += 1
    return count
