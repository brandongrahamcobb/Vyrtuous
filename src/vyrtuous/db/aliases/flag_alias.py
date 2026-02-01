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

    ARGS_MAP = {
        "alias_name": 1,
        "member": 2,
        "reason": 3
    }

    TABLE_NAME = "active_flags"

    @classmethod
    async def act_embed(cls, information, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been flagged",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Reason:** {infraction_information['infraction_reason']}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, information, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name}'s flag has been removed",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def enforce(
        cls, information, message, state
    ):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        flag = Flag(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
            reason=information["reason"],
        )
        await flag.create()
        bot = DiscordBot.get_instance()
        cog = bot.get_cog("ChannelEventListeners")
        cog.flags.append(flag)
        await Streaming.send_entry(
            alias=information['alias'],
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            duration=information["duration"],
            is_channel_scope=False,
            is_modification=information["modification"],
            member=member,
            message=message,
            reason=information["reason"],
        )
        embed = await FlagAlias.act_embed(
            information=information, source=message
        )
        return await state.end(success=embed)


    @classmethod
    async def undo(
        cls, information, message, state
    ):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        await Flag.delete(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
        )
        bot = DiscordBot.get_instance()
        cog = bot.get_cog("ChannelEventListeners")
        for flag in cog.flags:
            if flag.channel_snowflake == information["snowflake_kwargs"]["channel_snowflake"]:
                cog.flags.remove(flag)
                break
        await Streaming.send_entry(
            alias=information['alias'],
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            duration="",
            is_channel_scope=False,
            is_modification=True,
            member=member,
            message=message,
            reason="No reason provided.",
        )
        embed = await FlagAlias.undo_embed(
            information=information, source=message
        )
        return await state.end(success=embed)
