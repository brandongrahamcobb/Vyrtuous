"""!/bin/python3
clear_service.py The purpose of this program is to service the clear command.

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

import discord

from vyrtuous.aliases import (
    alias_service,
    unban_alias_service,
    unflag_alias_service,
    untext_mute_alias_service,
    unvoice_mute_alias_service,
)
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import PermissionState
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.channels import automute_channel_service, video_channel_service
from vyrtuous.utils.moderation import (
    ban_service,
    flag_service,
    text_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.permissions import permission_service
from vyrtuous.utils.tracking import stream_service

INFRACTION_MODELS = [
    ban_service.MODEL,
    flag_service.MODEL,
    text_mute_service.MODEL,
    voice_mute_service.MODEL,
]
CHANNEL_MODELS = [
    automute_channel_service.MODEL,
    stream_service.MODEL,
    video_channel_service.MODEL,
]
ALIAS_MODEL = [alias_service.MODEL]


async def clear(
    author_snowflake: int,
    category: str,
    guild_snowflake: int,
    message_snowflake: int,
    message_channel_snowflake: int,
    obj,
    view,
    *,
    target: str = "click",
):
    bot: DiscordBot = DiscordBot.get_instance()
    permission_state: PermissionState = bot.registry.get(PermissionState)
    if isinstance(obj, discord.Member):
        COMBINED_LIST = INFRACTION_MODELS
        await permission_service.has_equal_or_lower_role(
            permission_state=permission_state,
            channel_snowflake=message_channel_snowflake,
            guild_snowflake=guild_snowflake,
            author_snowflake=author_snowflake,
            member_snowflake=obj.id,
        )
        await permission_service.any_group_has_permissions(
            permission_state=permission_state,
            member_snowflake=author_snowflake,
            requested=["command.clear.scope.member"],
        )
        msg = f"Deleted all associated {category} records for {obj.mention}."
        if view.result:
            for model in COMBINED_LIST:
                database_factory = DatabaseFactory(model)
                objects = await database_factory.select(
                    guild_snowflake=guild_snowflake,
                    member_snowflake=obj.id,
                    singular=False,
                )
                for value in objects:
                    if category == "all":
                        msg = f"Deleted all database information for {obj.mention}."
                        await database_factory.delete_by_cls(value)
                    match value.identifier:
                        case "ban":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=message_channel_snowflake,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.ban"],
                            )
                            if category == "ban":
                                await database_factory.delete_by_cls(value)
                            await unban_alias_service.unban(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=obj.id,
                            )
                            await unban_alias_service.log_unban(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=obj.id,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "flag":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=message_channel_snowflake,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.flag"],
                            )
                            if category == "flag":
                                await database_factory.delete_by_cls(value)
                            await unflag_alias_service.unflag(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=obj.id,
                            )
                            await unflag_alias_service.log_unflag(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=obj.id,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "tmute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=message_channel_snowflake,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.text-mute"],
                            )
                            if category == "tmute":
                                await database_factory.delete_by_cls(value)
                            await untext_mute_alias_service.untext_mute(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=obj.id,
                            )
                            await untext_mute_alias_service.log_untext_mute(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=obj.id,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "vmute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=message_channel_snowflake,
                                guild_snowflake=guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.voice-mute"],
                            )
                            if category == "vmute":
                                await database_factory.delete_by_cls(value)
                            match target:
                                case "all":
                                    await permission_service.has_permissions(
                                        permission_state=permission_state,
                                        channel_snowflake=message_channel_snowflake,
                                        guild_snowflake=guild_snowflake,
                                        member_snowflake=author_snowflake,
                                        requested=["command.clear.target.all"],
                                    )
                                    targets = ["auto", "click", "command", "server"]
                                    for target in targets:
                                        await unvoice_mute_alias_service.unvoice_mute(
                                            channel_snowflake=value.channel_snowflake,
                                            guild_snowflake=guild_snowflake,
                                            member_snowflake=obj.id,
                                            target=target,
                                        )
                                        await unvoice_mute_alias_service.log_unvoice_mute(
                                            author_snowflake=author_snowflake,
                                            channel_snowflake=value.channel_snowflake,
                                            display=False,
                                            guild_snowflake=guild_snowflake,
                                            is_channel_scope=False,
                                            member_snowflake=obj.id,
                                            message_snowflake=message_snowflake,
                                            message_channel_snowflake=message_channel_snowflake,
                                            reason="Clear command.",
                                            target=target,
                                        )
                                case _:
                                    await permission_service.has_permissions(
                                        permission_state=permission_state,
                                        channel_snowflake=message_channel_snowflake,
                                        guild_snowflake=guild_snowflake,
                                        member_snowflake=author_snowflake,
                                        requested=[f"command.clear.target.{target}"],
                                    )
                                    await unvoice_mute_alias_service.unvoice_mute(
                                        channel_snowflake=value.channel_snowflake,
                                        guild_snowflake=guild_snowflake,
                                        member_snowflake=obj.id,
                                        target=target,
                                    )
                                    await unvoice_mute_alias_service.log_unvoice_mute(
                                        author_snowflake=author_snowflake,
                                        channel_snowflake=value.channel_snowflake,
                                        display=False,
                                        guild_snowflake=guild_snowflake,
                                        is_channel_scope=False,
                                        member_snowflake=obj.id,
                                        message_snowflake=message_snowflake,
                                        message_channel_snowflake=message_channel_snowflake,
                                        reason="Clear command.",
                                        target=target,
                                    )
    elif isinstance(obj, discord.abc.GuildChannel):
        COMBINED_LIST = INFRACTION_MODELS + CHANNEL_MODELS + ALIAS_MODEL
        await permission_service.any_group_has_permissions(
            permission_state=permission_state,
            member_snowflake=author_snowflake,
            requested=["command.clear.scope.channel"],
        )
        msg = f"Deleted all associated {category} records in {obj.mention}."
        if view.result:
            for model in COMBINED_LIST:
                database_factory = DatabaseFactory(model)
                objects = await database_factory.select(
                    channel_snowflake=obj.id,
                    guild_snowflake=obj.guild.id,
                    singular=False,
                )
                for value in objects:
                    await permission_service.has_equal_or_lower_role(
                        permission_state=permission_state,
                        channel_snowflake=obj.id,
                        guild_snowflake=value.guild_snowflake,
                        author_snowflake=author_snowflake,
                        member_snowflake=value.member_snowflake,
                    )
                    if category == "all":
                        msg = f"Deleted all database information for {obj.mention}."
                        await database_factory.delete_by_cls(value)
                    match value.identifier:
                        case "alias":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=obj.id,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.alias"],
                            )
                            if category == "alias":
                                await database_factory.delete_by_cls(value)
                            await alias_service.delete_alias(
                                alias_name=value.alias_name,
                                guild_snowflake=value.guild_snowflake,
                            )
                        case "automute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=obj.id,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.automute"],
                            )
                            if category == "automute":
                                await database_factory.delete_by_cls(value)
                            await automute_channel_service.toggle_automute(
                                author_snowflake=author_snowflake,
                                channel_snowflake=obj.id,
                                duration=None,
                                guild_snowflake=obj.guild.id,
                            )
                        case "ban":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=obj.id,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.ban"],
                            )
                            if category == "ban":
                                await database_factory.delete_by_cls(value)
                            await unban_alias_service.unban(
                                channel_snowflake=obj.id,
                                guild_snowflake=obj.guild.id,
                                member_snowflake=value.member_snowflake,
                            )
                            await unban_alias_service.log_unban(
                                author_snowflake=author_snowflake,
                                channel_snowflake=obj.id,
                                display=False,
                                guild_snowflake=obj.guild.id,
                                member_snowflake=value.member_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "flag":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=obj.id,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.flag"],
                            )
                            if category == "flag":
                                await database_factory.delete_by_cls(value)
                            await unflag_alias_service.unflag(
                                channel_snowflake=obj.id,
                                guild_snowflake=obj.guild.id,
                                member_snowflake=value.member_snowflake,
                            )
                            await unflag_alias_service.log_unflag(
                                author_snowflake=author_snowflake,
                                channel_snowflake=obj.id,
                                display=False,
                                guild_snowflake=obj.guild.id,
                                member_snowflake=value.member_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "stream":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=obj.id,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.stream"],
                            )
                            if category == "stream":
                                await database_factory.delete_by_cls(value)
                            bot: DiscordBot = DiscordBot.get_instance()
                            target_channel = bot.get_channel(value.channel_snowflake)
                            if isinstance(
                                target_channel,
                                (
                                    discord.VoiceChannel,
                                    discord.TextChannel,
                                    discord.StageChannel,
                                ),
                            ):
                                guild = bot.get_guild(value.guild_snowflake)
                                await stream_service.toggle_stream(
                                    target_channel=target_channel,
                                    source=guild,
                                )
                        case "tmute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=obj.id,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.text-mute"],
                            )
                            if category == "tmute":
                                await database_factory.delete_by_cls(value)
                            await untext_mute_alias_service.untext_mute(
                                channel_snowflake=obj.id,
                                guild_snowflake=obj.guild.id,
                                member_snowflake=value.member_snowflake,
                            )
                            await untext_mute_alias_service.log_untext_mute(
                                author_snowflake=author_snowflake,
                                channel_snowflake=obj.id,
                                display=False,
                                guild_snowflake=obj.guild.id,
                                member_snowflake=value.member_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "video":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=obj.id,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.video-channel"],
                            )
                            if category == "video":
                                await database_factory.delete_by_cls(value)
                            await video_channel_service.toggle_video_channel(
                                channel_snowflake=obj.id, guild_snowflake=obj.guild.id
                            )
                        case "vmute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=obj.id,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.voice_mute"],
                            )
                            if category == "vmute":
                                await database_factory.delete_by_cls(value)
                            match target:
                                case "all":
                                    await permission_service.has_permissions(
                                        permission_state=permission_state,
                                        channel_snowflake=obj.id,
                                        guild_snowflake=value.guild_snowflake,
                                        member_snowflake=author_snowflake,
                                        requested=[
                                            "command.clear.target.auto",
                                            "command.clear.target.click",
                                            "command.clear.target.command",
                                            "command.clear.target.server",
                                        ],
                                    )
                                    targets = ["auto", "click", "command", "server"]
                                    for target in targets:
                                        await unvoice_mute_alias_service.unvoice_mute(
                                            channel_snowflake=obj.id,
                                            guild_snowflake=obj.guild.id,
                                            member_snowflake=value.member_snowflake,
                                            target=target,
                                        )
                                        await unvoice_mute_alias_service.log_unvoice_mute(
                                            author_snowflake=author_snowflake,
                                            channel_snowflake=obj.id,
                                            display=False,
                                            guild_snowflake=obj.guild.id,
                                            is_channel_scope=False,
                                            member_snowflake=value.member_snowflake,
                                            message_snowflake=message_snowflake,
                                            message_channel_snowflake=message_channel_snowflake,
                                            reason="Clear command.",
                                            target=target,
                                        )
                                case _:
                                    await permission_service.has_permissions(
                                        permission_state=permission_state,
                                        channel_snowflake=obj.id,
                                        guild_snowflake=value.guild_snowflake,
                                        member_snowflake=author_snowflake,
                                        requested=[f"command.clear.target.{target}"],
                                    )
                                    await unvoice_mute_alias_service.unvoice_mute(
                                        channel_snowflake=obj.id,
                                        guild_snowflake=obj.guild.id,
                                        member_snowflake=value.member_snowflake,
                                        target=target,
                                    )
                                    await unvoice_mute_alias_service.log_unvoice_mute(
                                        author_snowflake=author_snowflake,
                                        channel_snowflake=obj.id,
                                        display=False,
                                        guild_snowflake=obj.guild.id,
                                        is_channel_scope=False,
                                        member_snowflake=value.member_snowflake,
                                        message_snowflake=message_snowflake,
                                        message_channel_snowflake=message_channel_snowflake,
                                        reason="Clear command.",
                                        target=target,
                                    )
    elif isinstance(obj, discord.Guild):
        COMBINED_LIST = INFRACTION_MODELS + CHANNEL_MODELS + ALIAS_MODEL
        await permission_service.any_group_has_permissions(
            permission_state=permission_state,
            member_snowflake=author_snowflake,
            requested=["command.clear.scope.channel"],
        )
        msg = f"Deleted all associated {category} records in {obj.name}."
        if view.result:
            for model in COMBINED_LIST:
                database_factory = DatabaseFactory(model)
                objects = await database_factory.select(
                    guild_snowflake=obj.id,
                    singular=False,
                )
                for value in objects:
                    await permission_service.has_equal_or_lower_role(
                        permission_state=permission_state,
                        channel_snowflake=message_channel_snowflake,
                        guild_snowflake=obj.id,
                        author_snowflake=author_snowflake,
                        member_snowflake=value.member_snowflake,
                    )
                    if category == "all":
                        await permission_service.has_permissions(
                            permission_state=permission_state,
                            channel_snowflake=obj.id,
                            guild_snowflake=value.guild_snowflake,
                            member_snowflake=author_snowflake,
                            requested=["command.clear.scope.all"],
                        )
                        msg = f"Deleted all database information for {obj.name}."
                        await database_factory.delete_by_cls(value)
                    match value.identifier:
                        case "alias":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=obj.id,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.alias"],
                            )
                            if category == "alias":
                                await database_factory.delete_by_cls(value)
                            await alias_service.delete_alias(
                                alias_name=value.alias_name,
                                guild_snowflake=obj.id,
                            )
                        case "automute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflkae,
                                guild_snowflake=obj.id,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.automute"],
                            )
                            if category == "automute":
                                await database_factory.delete_by_cls(value)
                            await automute_channel_service.toggle_automute(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                duration=None,
                                guild_snowflake=obj.id,
                            )
                        case "ban":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflkae,
                                guild_snowflake=obj.id,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.ban"],
                            )
                            if category == "ban":
                                await database_factory.delete_by_cls(value)
                            await unban_alias_service.unban(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=obj.id,
                                member_snowflake=value.member_snowflake,
                            )
                            await unban_alias_service.log_unban(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=obj.id,
                                member_snowflake=value.member_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "flag":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflkae,
                                guild_snowflake=obj.id,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.flag"],
                            )
                            if category == "flag":
                                await database_factory.delete_by_cls(value)
                            await unflag_alias_service.unflag(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=obj.id,
                                member_snowflake=value.member_snowflake,
                            )
                            await unflag_alias_service.log_unflag(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=obj.id,
                                member_snowflake=value.channel_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "stream":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflkae,
                                guild_snowflake=obj.id,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.stream"],
                            )
                            if category == "stream":
                                await database_factory.delete_by_cls(value)
                            bot: DiscordBot = DiscordBot.get_instance()
                            target_channel = bot.get_channel(value.channel_snowflake)
                            if isinstance(
                                target_channel,
                                (
                                    discord.VoiceChannel,
                                    discord.TextChannel,
                                    discord.StageChannel,
                                ),
                            ):
                                guild = bot.get_guild(value.guild_snowflake)
                                await stream_service.toggle_stream(
                                    target_channel=target_channel,
                                    source=guild,
                                )
                        case "tmute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflkae,
                                guild_snowflake=obj.id,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.text-mute"],
                            )
                            if category == "tmute":
                                await database_factory.delete_by_cls(value)
                            await untext_mute_alias_service.untext_mute(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=obj.id,
                                member_snowflake=value.member_snowflake,
                            )
                            await untext_mute_alias_service.log_untext_mute(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=obj.id,
                                member_snowflake=value.member_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "video":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflkae,
                                guild_snowflake=obj.id,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.video-channel"],
                            )
                            if category == "video":
                                await database_factory.delete_by_cls(value)
                            await video_channel_service.toggle_video_channel(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=obj.id,
                            )
                        case "vmute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflkae,
                                guild_snowflake=obj.id,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.voice-mute"],
                            )
                            if category == "vmute":
                                await database_factory.delete_by_cls(value)
                            match target:
                                case "all":
                                    await permission_service.has_permissions(
                                        permission_state=permission_state,
                                        channel_snowflake=value.channel_snowflkae,
                                        guild_snowflake=obj.id,
                                        member_snowflake=author_snowflake,
                                        requested=[
                                            "command.clear.target.auto",
                                            "command.clear.target.click",
                                            "command.clear.target.command",
                                            "command.clear.target.server",
                                        ],
                                    )
                                    targets = ["auto", "click", "command", "server"]
                                    for target in targets:
                                        await unvoice_mute_alias_service.unvoice_mute(
                                            channel_snowflake=value.channel_snowflake,
                                            guild_snowflake=obj.id,
                                            member_snowflake=value.member_snowflake,
                                            target=target,
                                        )
                                        await unvoice_mute_alias_service.log_unvoice_mute(
                                            author_snowflake=author_snowflake,
                                            channel_snowflake=value.channel_snowflake,
                                            display=False,
                                            guild_snowflake=obj.id,
                                            is_channel_scope=False,
                                            member_snowflake=value.member_snowflake,
                                            message_snowflake=message_snowflake,
                                            message_channel_snowflake=message_channel_snowflake,
                                            reason="Clear command.",
                                            target=target,
                                        )
                                case _:
                                    await permission_service.has_permissions(
                                        permission_state=permission_state,
                                        channel_snowflake=value.channel_snowflkae,
                                        guild_snowflake=obj.id,
                                        member_snowflake=author_snowflake,
                                        requested=[f"command.clear.target.{target}"],
                                    )
                                    await unvoice_mute_alias_service.unvoice_mute(
                                        channel_snowflake=value.channel_snowflake,
                                        guild_snowflake=obj.id,
                                        member_snowflake=value.member_snowflake,
                                        target=target,
                                    )
                                    await unvoice_mute_alias_service.log_unvoice_mute(
                                        author_snowflake=author_snowflake,
                                        channel_snowflake=value.channel_snowflake,
                                        display=False,
                                        guild_snowflake=obj.id,
                                        is_channel_scope=False,
                                        member_snowflake=value.member_snowflake,
                                        message_snowflake=message_snowflake,
                                        message_channel_snowflake=message_channel_snowflake,
                                        reason="Clear command.",
                                        target=target,
                                    )
    elif obj == "all" and permission_service.is_sysadmin(
        member_snowflake=author_snowflake
    ):
        COMBINED_LIST = INFRACTION_MODELS + CHANNEL_MODELS + ALIAS_MODEL
        msg = f"Deleted all associated {category} database entries."
        if view.result:
            for model in COMBINED_LIST:
                database_factory = DatabaseFactory(model)
                objects = await database_factory.select(
                    singular=False,
                )
                for value in objects:
                    if category == "all":
                        await permission_service.has_permissions(
                            permission_state=permission_state,
                            channel_snowflake=value.channel_snowflake,
                            guild_snowflake=value.guild_snowflake,
                            member_snowflake=author_snowflake,
                            requested=["command.clear.scope.all"],
                        )
                        msg = f"Deleted all database information."
                        await database_factory.delete_by_cls(value)
                    match value.identifier:
                        case "alias":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.alias"],
                            )
                            if category == "alias":
                                await database_factory.delete_by_cls(value)
                            await alias_service.delete_alias(
                                alias_name=value.alias_name,
                                guild_snowflake=value.guild_snowflake,
                            )
                        case "automute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.automute"],
                            )
                            if category == "automute":
                                await database_factory.delete_by_cls(value)
                            await automute_channel_service.toggle_automute(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                duration=None,
                                guild_snowflake=value.guild_snowflake,
                            )
                        case "ban":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.ban"],
                            )
                            if category == "ban":
                                await database_factory.delete_by_cls(value)
                            await unban_alias_service.unban(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=value.member_snowflake,
                            )
                            await unban_alias_service.log_unban(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=value.member_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "flag":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.flag"],
                            )
                            if category == "flag":
                                await database_factory.delete_by_cls(value)
                            await unflag_alias_service.unflag(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=value.member_snowflake,
                            )
                            await unflag_alias_service.log_unflag(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=value.member_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "stream":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.stream"],
                            )
                            if category == "stream":
                                await database_factory.delete_by_cls(value)
                            target_channel = bot.get_channel(value.channel_snowflake)
                            if isinstance(
                                target_channel,
                                (
                                    discord.VoiceChannel,
                                    discord.TextChannel,
                                    discord.StageChannel,
                                ),
                            ):
                                guild = bot.get_guild(value.guild_snowflake)
                                await stream_service.toggle_stream(
                                    target_channel=target_channel,
                                    source=guild,
                                )
                        case "tmute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.text-mute"],
                            )
                            if category == "tmute":
                                await database_factory.delete_by_cls(value)
                            await untext_mute_alias_service.untext_mute(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=value.member_snowflake,
                            )
                            await untext_mute_alias_service.log_untext_mute(
                                author_snowflake=author_snowflake,
                                channel_snowflake=value.channel_snowflake,
                                display=False,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=value.member_snowflake,
                                message_snowflake=message_snowflake,
                                message_channel_snowflake=message_channel_snowflake,
                                reason="Clear command.",
                            )
                        case "video":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.video-channel"],
                            )
                            if category == "video":
                                await database_factory.delete_by_cls(value)
                            await video_channel_service.toggle_video_channel(
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                            )
                        case "vmute":
                            await permission_service.has_permissions(
                                permission_state=permission_state,
                                channel_snowflake=value.channel_snowflake,
                                guild_snowflake=value.guild_snowflake,
                                member_snowflake=author_snowflake,
                                requested=["command.clear.category.voice-mute"],
                            )
                            if category == "vmute":
                                await database_factory.delete_by_cls(value)
                            match target:
                                case "all":
                                    await permission_service.has_permissions(
                                        permission_state=permission_state,
                                        channel_snowflake=value.channel_snowflake,
                                        guild_snowflake=value.guild_snowflake,
                                        member_snowflake=author_snowflake,
                                        requested=[
                                            "command.clear.target.auto",
                                            "command.clear.target.click",
                                            "command.clear.target.command",
                                            "command.clear.target.server",
                                        ],
                                    )
                                    targets = ["auto", "click", "command", "server"]
                                    for target in targets:
                                        await unvoice_mute_alias_service.unvoice_mute(
                                            channel_snowflake=value.channel_snowflake,
                                            guild_snowflake=value.guild_snowflake,
                                            member_snowflake=value.member_snowflake,
                                            target=target,
                                        )
                                        await unvoice_mute_alias_service.log_unvoice_mute(
                                            author_snowflake=author_snowflake,
                                            channel_snowflake=value.channel_snowflake,
                                            display=False,
                                            guild_snowflake=value.guild_snowflake,
                                            is_channel_scope=False,
                                            member_snowflake=value.member_snowflake,
                                            message_snowflake=message_snowflake,
                                            message_channel_snowflake=message_channel_snowflake,
                                            reason="Clear command.",
                                            target=target,
                                        )
                                case _:
                                    await permission_service.has_permissions(
                                        permission_state=permission_state,
                                        channel_snowflake=value.channel_snowflake,
                                        guild_snowflake=value.guild_snowflake,
                                        member_snowflake=author_snowflake,
                                        requested=[f"command.clear.target.{target}"],
                                    )
                                    await unvoice_mute_alias_service.unvoice_mute(
                                        channel_snowflake=value.channel_snowflake,
                                        guild_snowflake=value.guild_snowflake,
                                        member_snowflake=value.member_snowflake,
                                        target=target,
                                    )
                                    await unvoice_mute_alias_service.log_unvoice_mute(
                                        author_snowflake=author_snowflake,
                                        channel_snowflake=value.channel_snowflake,
                                        display=False,
                                        guild_snowflake=value.guild_snowflake,
                                        is_channel_scope=False,
                                        member_snowflake=value.member_snowflake,
                                        message_snowflake=message_snowflake,
                                        message_channel_snowflake=message_channel_snowflake,
                                        reason="Clear command.",
                                        target=target,
                                    )
    else:
        if obj is None:
            msg = f"Clear command cancelled."
        else:
            msg = f"Invalid target ({obj})."
    return msg
