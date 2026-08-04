"""!/bin/python3
snowflake_context.py  The purpose of this program is to provide a struct for snowflake context.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MemberState
from vyrtuous.utils.errors.error import ChannelNotFound, GuildNotFound, MemberNotFound


class SnowflakeContext:
    def __init__(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
    ):
        self.bot: DiscordBot = DiscordBot.get_instance()
        if guild_snowflake:
            self.guild_snowflake = guild_snowflake
            self.guild = self.bot.get_guild(guild_snowflake)
            if self.guild is None:
                raise GuildNotFound(str(guild_snowflake))
            if channel_snowflake:
                self.channel_snowflake = channel_snowflake
                self.channel = self.guild.get_channel(channel_snowflake)
                if self.channel is None:
                    raise ChannelNotFound(str(channel_snowflake))
            if member_snowflake:
                self.member_snowflake = member_snowflake
                member = self.guild.get_member(member_snowflake)
                if member is None:
                    simplified_member = self.bot.registry.get(MemberState).active.get(
                        member_snowflake, None
                    )
                    if simplified_member is None:
                        raise MemberNotFound(str(member_snowflake))
