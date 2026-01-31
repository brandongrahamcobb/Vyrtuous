"""flag.py The purpose of this program is to inherit from DatabaseFactory to provide the flag moderation.

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

from vyrtuous.db.infractions.flag import Flag
from vyrtuous.db.mgmt.stream import Streaming
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.emojis import get_random_emoji


class FlagAlias(Alias):

    identifier = "flag"

    ACT = "flag"
    UNDO = "unflag"

    @classmethod
    async def act_embed(cls, infraction_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(infraction_information["infraction_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(infraction_information["infraction_member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been flagged",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Reason:** {infraction_information['infraction_reason']}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, infraction_information, source, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(infraction_information["infraction_channel_snowflake"])
        author = resolve_author(source=source)
        member = source.guild.get_member(infraction_information["infraction_member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name}'s flag has been removed",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def enforce(
        cls, alias, infraction_information, member, message, state
    ):
        flag = Flag(
            channel_snowflake=infraction_information["infraction_channel_snowflake"],
            guild_snowflake=infraction_information["infraction_guild_snowflake"],
            member_snowflake=infraction_information["infraction_member_snowflake"],
            reason=infraction_information["infraction_reason"],
        )
        await flag.create()
        bot = DiscordBot.get_instance()
        cog = bot.get_cog("ChannelEventListeners")
        cog.flags.append(flag)
        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=infraction_information["infraction_channel_snowflake"],
            duration=infraction_information["infraction_duration"],
            is_channel_scope=False,
            is_modification=infraction_information["infraction_modification"],
            member=member,
            message=message,
            reason=infraction_information["infraction_reason"],
        )
        embed = await Flag.act_embed(
            infraction_information=infraction_information, source=message
        )
        return await state.end(success=embed)


    @classmethod
    async def undo(
        cls, alias, infraction_information, member, message, state
    ):
        await Flag.delete(
            channel_snowflake=infraction_information["infraction_channel_snowflake"],
            guild_snowflake=infraction_information["infraction_guild_snowflake"],
            member_snowflake=infraction_information["infraction_member_snowflake"],
        )
        bot = DiscordBot.get_instance()
        cog = bot.get_cog("ChannelEventListeners")
        for flag in cog.flags:
            if flag.channel_snowflake == infraction_information["infraction_channel_snowflake"]:
                cog.flags.remove(flag)
                break
        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=infraction_information["infraction_channel_snowflake"],
            duration="",
            is_channel_scope=False,
            is_modification=infraction_information["infraction_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )
        embed = await Flag.undo_embed(
            infraction_information=infraction_information, source=message
        )
        return await state.end(success=embed)
