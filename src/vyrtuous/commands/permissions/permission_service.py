from typing import Dict, Tuple, Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.commands.errors import HasEqualOrLowerRole
from vyrtuous.db.infractions.ban.ban import Ban
from vyrtuous.db.infractions.flag.flag import Flag
from vyrtuous.db.infractions.tmute.text_mute import TextMute
from vyrtuous.db.infractions.vmute.voice_mute import VoiceMute
from vyrtuous.db.roles.admin.administrator_service import (
    NotAdministrator,
    is_administrator,
    is_administrator_at_all
)
from vyrtuous.db.roles.coord.coordinator_service import (
    NotCoordinator,
    is_coordinator,
    is_coordinator_at_all,
)
from vyrtuous.db.roles.dev.developer_service import NotDeveloper, is_developer
from vyrtuous.db.roles.mod.moderator_service import (
    NotModerator,
    is_moderator,
    is_moderator_at_all,
)
from vyrtuous.db.roles.owner.guild_owner_service import NotGuildOwner, is_guild_owner, is_guild_owner_at_all
from vyrtuous.db.roles.sysadmin.sysadmin_service import NotSysadmin, is_sysadmin
from vyrtuous.inc.helpers import CHUNK_SIZE, PERMISSION_TYPES, TARGET_PERMISSIONS
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_channels,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.logger import logger


class PermissionService:

    invincible_members: Dict[Tuple[int, int], bool] = {}
    lines, pages = [], []
    state: bool = False

    @classmethod
    async def unrestrict(cls, guild_snowflake, member_snowflake):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        kwargs = {
            "guild_snowflake": int(guild_snowflake),
            "member_snowflake": member_snowflake,
        }
        bans = await Ban.select(**kwargs)
        text_mutes = await TextMute.select(**kwargs)
        voice_mutes = await VoiceMute.select(**kwargs)
        if bans:
            for ban in bans:
                channel = guild.get_channel(ban.channel_snowflake)
                if channel:
                    try:
                        await channel.set_permissions(member, overwrite=None)
                    except discord.Forbidden:
                        logger.warning(
                            f"Unable to unban {member.name} ({member.id}) in {channel.name} ({channel.id})."
                        )
        if text_mutes:
            for text_mute in text_mutes:
                channel = guild.get_channel(text_mute.channel_snowflake)
                if channel:
                    try:
                        await channel.set_permissions(member, send_messages=True)
                    except discord.Forbidden:
                        logger.warning(
                            f"Unable to untmute {member.name} ({member.id}) in {channel.name} ({channel.id})."
                        )
        if voice_mutes:
            for voice_mute in voice_mutes:
                channel = guild.get_channel(voice_mute.channel_snowflake)
                if channel and member.voice and member.voice.mute:
                    await member.edit(mute=False)
        await Ban.delete(**kwargs)
        await Flag.delete(**kwargs)
        await TextMute.delete(**kwargs)
        await VoiceMute.delete(**kwargs)

    @classmethod
    def add_invincible_member(cls, guild_snowflake: int, member_snowflake: int):
        cls.invincible_members[(guild_snowflake, member_snowflake)] = True

    @classmethod
    def get_invincible_members(cls):
        return cls.invincible_members

    @classmethod
    def remove_invincible_member(cls, guild_snowflake: int, member_snowflake: int):
        cls.invincible_members.pop((guild_snowflake, member_snowflake), None)

    @classmethod
    def toggle_enabled(cls):
        cls.state = not cls.state
        return cls.state

    @classmethod
    async def build_clean_dictionary(cls, channel_objs, is_at_home, me):
        dictionary = {}
        for channel in channel_objs:
            permissions = channel.permissions_for(me)
            missing = []
            for permission in TARGET_PERMISSIONS:
                if not getattr(permissions, permission):
                    missing.append(permission)
            if not missing:
                continue
            dictionary.setdefault(channel.guild.id, {"channels": {}})
            dictionary[channel.guild.id]["channels"].setdefault(channel.id, {})
            dictionary[channel.guild.id]["channels"][channel.id].update(
                {"permissions": missing}
            )
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_channels:
                PermissionService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
            if skipped_guilds:
                PermissionService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, is_at_home, channel_objs, default_kwargs):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        channel_snowflake = default_kwargs.get("channel_snowflake", None)
        guild_snowflake = default_kwargs.get("guild_snowflake", None)
        guild = bot.get_guild(guild_snowflake)
        title = f"{get_random_emoji()} {bot.user.display_name} Missing Permissions"

        dictionary = await PermissionService.build_clean_dictionary(
            channel_objs=channel_objs, is_at_home=is_at_home, me=guild.me
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                PermissionService.lines.append(f"Channel: {channel.mention}")
                for section_name, permissions in channel_data.items():
                    for permission in permissions:
                        PermissionService.lines.append(f"  ↳ {permission}")
                field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(PermissionService.lines),
                        inline=False,
                    )
                    embed = flush_page(
                        embed, PermissionService.pages, title, guild.name
                    )
                    PermissionService.lines = []
            if PermissionService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(PermissionService.lines),
                    inline=False,
                )
            PermissionService.pages.append(embed)
        return PermissionService.pages

    @classmethod
    async def check(
        cls,
        updated_kwargs,
        lowest_role: str,
    ) -> str:
        verifications = (
            ("Sysadmin", is_sysadmin),
            ("Developer", is_developer),
            ("Guild Owner", is_guild_owner),
            ("Administrator", is_administrator),
            ("Coordinator", is_coordinator),
            ("Moderator", is_moderator),
        )
        channel_snowflake = updated_kwargs.get("channel_snowflake", None)
        guild_snowflake = updated_kwargs.get("guild_snowflake", None)
        member_snowflake = updated_kwargs.get("member_snowflake", None)
        passed_lowest = False
        for role_name, verify in verifications:
            try:
                if role_name in ("Sysadmin", "Developer"):
                    if await verify(member_snowflake=int(member_snowflake)):
                        return role_name
                elif role_name in ("Guild Owner", "Administrator"):
                    if await verify(
                        guild_snowflake=int(guild_snowflake),
                        member_snowflake=int(member_snowflake),
                    ):
                        return role_name
                else:
                    if await verify(
                        channel_snowflake=int(channel_snowflake),
                        guild_snowflake=int(guild_snowflake),
                        member_snowflake=int(member_snowflake),
                    ):
                        return role_name
            except commands.CheckFailure:
                if lowest_role is not None and passed_lowest:
                    raise
            if role_name == lowest_role:
                passed_lowest = True
        return "Everyone"

    @classmethod
    async def has_equal_or_lower_role(
        cls,
        updated_kwargs,
        member_snowflake: int,
    ) -> bool:
        sender_name = await PermissionService.resolve_highest_role(**updated_kwargs)
        sender_rank = PERMISSION_TYPES.index(sender_name)
        new_kwargs = updated_kwargs.copy()
        new_kwargs.update({"member_snowflake": member_snowflake})
        target_name = await PermissionService.resolve_highest_role(**new_kwargs)
        target_rank = PERMISSION_TYPES.index(target_name)
        if sender_rank <= target_rank:
            raise HasEqualOrLowerRole(PERMISSION_TYPES[target_rank])
        return sender_name

    @classmethod
    async def resolve_highest_role(
        cls,
        channel_snowflake: int,
        member_snowflake: int,
        guild_snowflake: int,
    ):
        try:
            if await is_sysadmin(member_snowflake=int(member_snowflake)):
                return "Sysadmin"
        except NotSysadmin as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_developer(member_snowflake=int(member_snowflake)):
                return "Developer"
        except NotDeveloper as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_guild_owner(
                guild_snowflake=int(guild_snowflake),
                member_snowflake=int(member_snowflake),
            ):
                return "Guild Owner"
        except NotGuildOwner as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_administrator(
                guild_snowflake=int(guild_snowflake),
                member_snowflake=int(member_snowflake),
            ):
                return "Administrator"
        except NotAdministrator as e:
            logger.warning(str(e).capitalize())
        if channel_snowflake:
            try:
                if await is_coordinator(
                    channel_snowflake=int(channel_snowflake),
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=int(member_snowflake),
                ):
                    return "Coordinator"
            except NotCoordinator as e:
                logger.warning(str(e).capitalize())
            try:
                if await is_moderator(
                    channel_snowflake=int(channel_snowflake),
                    guild_snowflake=int(guild_snowflake),
                    member_snowflake=int(member_snowflake),
                ):
                    return "Moderator"
            except NotModerator as e:
                logger.warning(str(e).capitalize())
        return "Everyone"

    @classmethod
    async def resolve_highest_role_at_all(
        cls,
        member_snowflake: int,
    ):
        try:
            if await is_sysadmin(member_snowflake=int(member_snowflake)):
                return "Sysadmin"
        except NotSysadmin as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_developer(member_snowflake=int(member_snowflake)):
                return "Developer"
        except NotDeveloper as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_guild_owner_at_all(
                member_snowflake=int(member_snowflake),
            ):
                return "Guild Owner"
        except NotGuildOwner as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_administrator_at_all(
                member_snowflake=int(member_snowflake),
            ):
                return "Administrator"
        except NotAdministrator as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_coordinator_at_all(
                member_snowflake=int(member_snowflake),
            ):
                return "Coordinator"
        except NotCoordinator as e:
            logger.warning(str(e).capitalize())
        try:
            if await is_moderator_at_all(
                member_snowflake=int(member_snowflake),
            ):
                return "Moderator"
        except NotModerator as e:
            logger.warning(str(e).capitalize())
        return "Everyone"
