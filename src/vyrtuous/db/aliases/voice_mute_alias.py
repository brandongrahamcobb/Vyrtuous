"""voice_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the voice mute moderation.

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

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.infractions.voice_mute import VoiceMute
from vyrtuous.db.mgmt.stream import Streaming
from vyrtuous.utils.author import resolve_author
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class VoiceMuteAlias(Alias):
    
    identifier = "vmute"

    ACT = "vmute"
    UNDO = "unvmute"

    ARGS_MAP = {
        "alias_name": 1,
        "member": 2,
        "duration": 3,
        "reason": 4
    }
    
    TABLE_NAME = "active_voice_mutes"

    @classmethod
    async def act_embed(cls, information, **kwargs):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(information["snowflake_kwargs"]["channel_snowflake"])
        guild = bot.get_guild(information["snowflake_kwargs"]["guild_snowflake"])
        member = guild.get_member(information["snowflake_kwargs"]["member_snowflake"])
        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name} has been voice-muted",
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
            title=f"{get_random_emoji()} "
            f"{member.display_name}'s voice-mute has been removed",
            description=(
                f"**By:** {author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
            ),
            color=discord.Color.blue(),
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
        voice_mute = VoiceMute(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            expires_in=information["expires_in"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
            reason=information["reason"],
            target="user",
        )
        await voice_mute.create()
        is_channel_scope = False
        if member.voice and member.voice.channel:
            if (
                member.voice.channel.id
                == information["snowflake_kwargs"]["channel_snowflake"]
            ):
                is_channel_scope = True
                try:
                    await member.edit(
                        mute=True, reason=information["reason"]
                    )
                except discord.Forbidden as e:
                    return await state.end(error=str(e).capitalize())
        await Streaming.send_entry(
            alias=information['alias'],
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            duration=information["duration"],
            is_channel_scope=is_channel_scope,
            member=member,
            message=message,
            reason=information["reason"],
        )
        embed = await VoiceMuteAlias.act_embed(
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
        await VoiceMute.delete(
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            guild_snowflake=information["snowflake_kwargs"]["guild_snowflake"],
            member_snowflake=information["snowflake_kwargs"]["member_snowflake"],
        )
        is_channel_scope = False
        if member.voice and member.voice.channel:
            try:
                is_channel_scope = True
                await member.edit(mute=False)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())
        await Streaming.send_entry(
            alias=information['alias'],
            channel_snowflake=information["snowflake_kwargs"]["channel_snowflake"],
            is_channel_scope=is_channel_scope,
            is_modification=True,
            member=member,
            message=message,
        )
        embed = await VoiceMuteAlias.undo_embed(
            information=information, source=message
        )
        return await state.end(success=embed)
