"""!/bin/python3
channel_event_listeners.py A discord.py cog containing channel event listeners for the Vyrtuous bot.

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
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import ChannelState, MemberState
from vyrtuous.utils.moderation import (ban_service, flag_service,
                                       server_mute_service, voice_mute_service)
from vyrtuous.utils.rooms import (automute_room_service,
                                  temporary_room_service, video_room_service)
from vyrtuous.utils.users import moderator_service


class ChannelEventListeners(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    @commands.Cog.listener()
    async def on_guild_channel_grant(self, channel: discord.abc.GuildChannel):
        guild = channel.guild
        name = channel.name
        for c in guild.channels:
            if c.id != channel.id and c.name == name:
                return
        room = self.__bot.registry.get(ChannelState).deleted.pop(name, None)
        await temporary_room_service.migrate_temporary_room(
            channel=channel, old_name=room.name
        )

    @commands.Cog.listener()
    async def on_guild_channel_delete(self, channel: discord.abc.GuildChannel):
        await temporary_room_service.add_deleted_room(channel=channel)

    @commands.Cog.listener()
    async def on_guild_channel_update(self, before, after):
        if before.name == after.name:
            return
        await temporary_room_service.rename_room(before=before, after=after)

    @commands.Cog.listener()
    async def on_voice_state_update(
        self,
        member: discord.Member,
        before: discord.VoiceState,
        after: discord.VoiceState,
    ):
        if member.bot:
            return
        self.__bot.registry.get(MemberState).active[member.id] = (
            member.display_name,
            datetime.now(timezone.utc),
        )
        if not after.channel:
            return
        if before.channel == after.channel:
            if before.mute == after.mute:
                if before.self_mute == after.self_mute:
                    return
                return
        await ban_service.is_banned_then_kick_and_reset_cooldown(
            channel=after.channel, member=member
        )
        if await server_mute_service.is_server_muted(
            channel=after.channel, member=member
        ):
            return
        if video_room_service.is_active_video_room(channel=after.channel):
            await video_room_service.update_video_room_tasks(
                after=after, before=before, member=member
            )
        duration_value = "1h"
        if (
            after.channel.guild.id,
            member.id,
        ) in self.__bot.registry.get(MemberState).invincible:
            # embed = discord.Embed(
            #     title=f"{member.display_name} is a hero!",
            #     description=f"{member.display_name} cannot be muted.",
            #     color=discord.Color.gold(),
            # )
            await voice_mute_service.unmute(
                channel=after.channel, member=member, target="user"
            )
            return
            # embed.set_thumbnail(url=member.display_avatar.url)
            # return await after.channel.send(embed=embed)
        elif (
            await automute_room_service.is_active_stage_room(channel=after.channel)
            and await moderator_service.resolve_highest_role(
                channel_snowflake=after.channel.id,
                member_snowflake=member.id,
                guild_snowflake=after.channel.guild.id,
            )
            in "Everyone"
        ):
            target = "room"
            await voice_mute_service.mute(
                channel=after.channel,
                duration_value=duration_value,
                member=member,
                target=target,
            )
        elif before.channel != after.channel:
            if await voice_mute_service.is_voice_muted(
                channel=after.channel, member=member
            ):
                target = "user"
                await voice_mute_service.mute(
                    channel=after.channel,
                    duration_value=duration_value,
                    member=member,
                    target=target,
                )
            elif after.mute:
                target = "user"
                await voice_mute_service.unmute(
                    channel=after.channel, member=member, target=target
                )
            await flag_service.warn(channel=after.channel, member=member)
        elif before.channel == after.channel:
            target = "user"
            if before.mute and not after.mute:
                await voice_mute_service.delete(
                    channel_snowflake=after.channel.id,
                    guild_snowflake=after.channel.guild.id,
                    member_snowflake=member.id,
                )
            else:
                await voice_mute_service.mute(
                    channel=after.channel,
                    duration_value=duration_value,
                    member=member,
                    target=target,
                )


async def setup(bot: DiscordBot):
    await bot.add_cog(ChannelEventListeners(bot))
