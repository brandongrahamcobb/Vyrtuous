"""text_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the text mute moderation.

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

from datetime import datetime, timezone


import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.infractions.text_mute import TextMute
from vyrtuous.db.mgmt.stream import Streaming
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class TextMuteAlias(Alias):

    identifier = "tmute"

    ACT = "tmute"
    UNDO = "untmute"

    ARGS_MAP = {
        "alias_name": 1,
        "member": 2,
        "duration": 3,
        "reason": 4
    }

    TABLE_NAME = "active_text_mutes"

    @classmethod
    async def act_embed(cls, information, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Text-Muted",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {information['infraction_duration']}\n"
                f"**Reason:** {information['infraction_reason']}"
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
            title=f"{get_random_emoji()} " f"{member.display_name} has been Unmuted",
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
        text_mute = TextMute(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            expires_in=information["expires_in"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
            role_snowflake=information["snowflake_kwargs"]["role_snowflake"],
            reason=information["reason"],
        )
        await text_mute.create()
        channel = message.guild.get_channel(
            information["snowflake_kwargs"]["channel_snowflake"]
        )
        if channel:
            try:
                await channel.set_permissions(
                    target=member,
                    send_messages=False,
                    add_reactions=False,
                    reason=information["reason"],
                )
            except discord.Forbidden as e:
                logger.error(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        await Streaming.send_entry(
            alias=information['alias'],
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            duration=information["duration"],
            member=member,
            message=message,
            reason=information["reason"],
        )
        embed = await TextMuteAlias.act_embed(
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
        await TextMute.delete(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
        )
        channel = message.guild.get_channel(
            information["snowflake_kwargs"]["channel_snowflake"]
        )
        if channel:
            try:
                await channel.set_permissions(
                    target=member,
                    send_messages=None,
                    add_reactions=None,
                    reason=information["reason"],
                )
            except discord.Forbidden as e:
                logger.error(str(e).capitalize())
                return await state.end(error=str(e).capitalize())
        await Streaming.send_entry(
            alias=information['alias'],
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            is_modification=True,
            member=member,
            message=message,
        )
        embed = await TextMuteAlias.undo_embed(
            information=information, source=message
        )
        return await state.end(success=embed)
