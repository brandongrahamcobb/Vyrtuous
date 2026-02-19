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

import time
from collections import defaultdict
from datetime import datetime, timedelta, timezone

import discord
from discord.ext import commands

from vyrtuous.alias.alias import Alias
from vyrtuous.alias.alias_service import AliasService
from vyrtuous.ban.ban import Ban
from vyrtuous.ban.ban_service import BanService
from vyrtuous.base.database_factory import DatabaseFactory
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cap.cap import Cap
from vyrtuous.cap.cap_service import CapService
from vyrtuous.coordinator.coordinator import Coordinator
from vyrtuous.coordinator.coordinator_service import CoordinatorService
from vyrtuous.duration.duration_service import DurationService
from vyrtuous.field.duration import DurationObject
from vyrtuous.flag.flag import Flag
from vyrtuous.flag.flag_service import FlagService
from vyrtuous.moderator.moderator import Moderator
from vyrtuous.moderator.moderator_service import ModeratorService
from vyrtuous.server_mute.server_mute import ServerMute
from vyrtuous.server_mute.server_mute_service import ServerMuteService
from vyrtuous.stage_room.stage import Stage
from vyrtuous.stage_room.stage_service import StageService
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.temporary_room.temporary_room import TemporaryRoom
from vyrtuous.temporary_room.temporary_room_service import TemporaryRoomService
from vyrtuous.text_mute.text_mute import TextMute
from vyrtuous.text_mute.text_mute_service import TextMuteService
from vyrtuous.utils.author_service import AuthorService
from vyrtuous.utils.dictionary_service import DictionaryService
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.hero_service import HeroService
from vyrtuous.utils.logger import logger
from vyrtuous.utils.message_service import PaginatorService
from vyrtuous.utils.permission_service import PermissionService
from vyrtuous.vegan.vegan_service import VeganService
from vyrtuous.video_room.video_room import VideoRoom
from vyrtuous.video_room.video_room_service import VideoRoomService
from vyrtuous.voice_mute.voice_mute import VoiceMute
from vyrtuous.voice_mute.voice_mute_service import VoiceMuteService


class ChannelEventListeners(commands.Cog):
    def __init__(self, bot: DiscordBot):
        self.__author_service = AuthorService()
        self.__bot = bot
        self.__database_factory = DatabaseFactory(bot=self.__bot)
        self.__dictionary_service = DictionaryService(bot=self.__bot)
        self.__emoji = Emojis()
        self.__duration_service = DurationService()
        self.__paginator_service = PaginatorService(bot=self.__bot)
        self.__moderator_service = ModeratorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__video_room_service = VideoRoomService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__stream_service = StreamService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            moderator_service=self.__moderator_service,
            paginator_service=self.__paginator_service,
        )
        self.__flag_service = FlagService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__voice_mute_service = VoiceMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
            stream_service=self.__stream_service,
        )
        self.__ban_service = BanService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__vegan_service = VeganService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__text_mute_service = TextMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
            stream_service=self.__stream_service,
        )
        self.__coordinator_service = CoordinatorService(
            author_service=self.__author_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__alias_service = AliasService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__stage_service = StageService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            moderator_service=self.__moderator_service,
        )
        self.__cap_service = CapService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            duration_service=self.__duration_service,
            emoji=self.__emoji,
        )
        self.__server_mute_service = ServerMuteService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
        )
        self.__hero_service = HeroService(
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            voice_mute_service=self.__voice_mute_service,
            ban_service=self.__ban_service,
            flag_service=self.__flag_service,
            text_mute_service=self.__text_mute_service,
        )
        self.__temporary_room_service = TemporaryRoomService(
            alias_service=self.__alias_service,
            bot=self.__bot,
            database_factory=self.__database_factory,
            dictionary_service=self.__dictionary_service,
            emoji=self.__emoji,
            cap_service=self.__cap_service,
            moderator_service=self.__moderator_service,
            stage_service=self.__stage_service,
            coordinator_service=self.__coordinator_service,
            voice_mute_service=self.__voice_mute_service,
            ban_service=self.__ban_service,
            flag_service=self.__flag_service,
            vegan_service=self.__vegan_service,
            text_mute_service=self.__text_mute_service,
        )

    async def cog_load(self):
        await self.__video_room_service.load_video_rooms_into_memory()
        await self.__flag_service.load_flags_into_memory()

    @commands.Cog.listener()
    async def on_guild_channel_grant(self, channel: discord.abc.GuildChannel):
        guild = channel.guild
        name = channel.name
        for c in guild.channels:
            if c.id != channel.id and c.name == name:
                return
        room = self.__temporary_room_service.deleted_rooms.pop(name, None)
        channel_dict = {
            "columns": {"channel_snowflake": channel.id, "guild_snowflake": guild.id},
            "id": channel.id,
            "name": name,
        }
        default_kwargs = {"guild_snowflake": guild.id}
        await self.__temporary_room_service.migrate_temporary_room(
            channel_dict=channel_dict, default_kwargs=default_kwargs, old_name=room.name
        )

    @commands.Cog.listener()
    async def on_guild_channel_delete(self, channel: discord.abc.GuildChannel):
        await self.__temporary_room_service.add_deleted_room(channel=channel)

    @commands.Cog.listener()
    async def on_guild_channel_update(self, before, after):
        if before.name == after.name:
            return
        await self.__temporary_room_service.rename_room(before=before, after=after)

    @commands.Cog.listener()
    async def on_voice_state_update(
        self,
        member: discord.Member,
        before: discord.VoiceState,
        after: discord.VoiceState,
    ):
        if member.bot:
            return
        if before.channel == after.channel:
            if before.mute == after.mute:
                if before.self_mute == after.self_mute:
                    return
                return
        await self.__ban_service.is_banned_then_kick_and_reset(
            channel=after.channel, member=member
        )
        if self.__server_mute_service.is_server_muted(
            channel=after.channel, member=member
        ):
            return
        await self.__video_room_service.update_video_room_tasks(
            after=after, before=before, member=member
        )
        expires_in = datetime.now(timezone.utc) + timedelta(hours=1)
        duration = self.__duration_service.from_expires_in(expires_in)
        if member.id in self.__hero_service.invincible_members:
            embed = discord.Embed(
                title=f"\u1f4aB {member.display_name} is a hero!",
                description=f"{member.display_name} cannot be muted.",
                color=discord.Color.gold(),
            )
            embed.set_thumbnail(url=member.display_avatar.url)
            return await after.channel.send(embed=embed)
        elif (
            await self.__stage_service.is_active_stage_room(channel=after.channel)
            and await self.__moderator_service.resolve_highest_role(
                channel_snowflake=after.channel.id,
                member_snowflake=member.id,
                guild_snowflake=after.channel.guild.id,
            )
            in "Everyone"
        ):
            target = "room"
            await self.__voice_mute_service.mute(
                channel=after.channel, duration=duration, member=member, target=target
            )
        elif await self.__voice_mute_service.is_voice_muted(
            channel=after.channel, member=member
        ):
            target = "user"
            if before.mute and not after.mute and before.channel == after.channel:
                await self.__voice_mute_service.unmute(
                    channel=after.channel, member=member, target=target
                )
            else:
                await self.__voice_mute_service.mute(
                    channel=after.channel,
                    duration=duration,
                    member=member,
                    target=target,
                )
        return await self.__flag_service.warn(channel=after.channel, member=member)


async def setup(bot: DiscordBot):
    await bot.add_cog(ChannelEventListeners(bot))
