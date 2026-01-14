"""aliases.py A discord.py cog containing command aliases for the Vyrtuous bot.

Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.actions.ban import Ban
from vyrtuous.database.actions.flag import Flag
from vyrtuous.database.actions.text_mute import TextMute
from vyrtuous.database.actions.voice_mute import VoiceMute
from vyrtuous.database.logs.history import History
from vyrtuous.database.roles.vegan import Vegan
from vyrtuous.database.settings.cap import Cap
from vyrtuous.properties.duration import DurationObject
from vyrtuous.service.logging_service import logger
from vyrtuous.service.resolution.role_service import resolve_role
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.invincibility import Invincibility


class Aliases(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.alias_help = {
            "ban": [
                "**member** (Required): Tag a member or include their ID",
                "**duration** (Optional): (+|-)duration(m|h|d)\n"
                "0 = permanent / 24h = default\n`+` to append, "
                "`-` to delete, `=` to overwrite reason",
                "**reason** (Optional): Reason (required for 7 days or more)",
            ],
            "vegan": ["**member** (Required): Tag a member or include their ID"],
            "carnist": ["**member** (Required): Tag a member or include their ID"],
            "unban": ["**member** (Required): Tag a member or include their ID"],
            "flag": [
                "**member** (Required): Tag a member or include their ID",
                "**reason** (Optional): Reason for flagging the user",
            ],
            "unflag": ["**member** (Required): Tag a member or include their ID"],
            "voice_mute": [
                "**member** (Required): Tag a member or include their ID",
                "**duration** (Optional): (+|-)duration(m|h|d)\n"
                "0 = permanent / 24h = default\n`+` to append, "
                "`-` to delete, `=` to overwrite reason",
                "**reason** (Optional): Reason (required for 7 days or more)",
            ],
            "unvoice_mute": ["**member** (Required): Tag a member or include their ID"],
            "text_mute": [
                "**member** (Required): Tag a member or include their ID",
                "**duration** (Required): (+|-)duration(m|h|d)\n"
                "0 = permanent / 24h = default\n`+` to append, "
                "`-` to delete, `=` to overwrite reason",
                "**reason** (Required): Reason (required for 7 days or more)",
            ],
            "untext_mute": ["**member** (Required): Tag a member or include their ID"],
            "role": [
                "**member** (Required): Tag a member or include their ID",
                "**role** (Required): Role to assign",
            ],
            "unrole": [
                "**member** (Required): Tag a member or include their ID",
                "**role** (Required): Role to remove",
            ],
        }
        self.alias_type_to_description = {
            "ban": "Bans a user from the server.",
            "vegan": "Verifies a user as going vegan.",
            "carnist": "Unverifies a user as going vegan.",
            "unban": "Unbans a user from the server.",
            "flag": "Flags a user for moderation review.",
            "unflag": "Removes a flag from a user.",
            "voice_mute": "Mutes a user in voice channels.",
            "unvoice_mute": "Unmutes a user in voice channels.",
            "text_mute": "Mutes a user in text channels.",
            "untext_mute": "Unmutes a user in text channels.",
            "role": "Assigns a role to a user.",
            "unrole": "Removes a role from a user.",
        }
        self.alias_type_to_permission_level = {
            "ban": "Moderator",
            "vegan": "Moderator",
            "carnist": "Moderator",
            "unban": "Moderator",
            "voice_mute": "Moderator",
            "unvoice_mute": "Moderator",
            "text_mute": "Moderator",
            "untext_mute": "Moderator",
            "flag": "Moderator",
            "unflag": "Moderator",
            "role": "Coordinator",
            "unrole": "Coordinator",
        }
        self.bot = bot
        self.invincible_members = Invincibility.get_invincible_members()

    async def handle_ban_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        try:
            bans= None
            is_channel_scope = False
            is_modification = False
            override = False

            reason = "No reason provided."
            if len(args) > 2:
                reason = " ".join(args[2:])

            cap_duration = generate_cap_duration(channel_snowflake=channel_obj.id, guild_snowflake=message.guild.id, moderation_type="ban")

            try:
                await channel_obj.set_permissions(
                    member_obj, view_channel=False, reason=reason
                )
            except discord.Forbidden as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
            if (
                member_obj.voice
                and member_obj.voice.channel
                and member_obj.voice.channel.id == channel_obj.id
            ):
                is_channel_scope = True
                try:
                    await member_obj.move_to(None, reason=reason)
                except discord.Forbidden as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")

            await History.send_entry(
                alias=alias,
                channel=channel_obj,
                duration=duration,
                is_channel_scope=is_channel_scope,
                is_modification=is_modification,
                member=member_obj,
                message=message,
                reason=reason,
            )

            embed = discord.Embed(
                title=f"{get_random_emoji()} "
                f"{member_obj.display_name} has been Banned",
                description=(
                    f"**By:** {message.author.mention}\n"
                    f"**User:** {member_obj.mention}\n"
                    f"**Channel:** {channel_obj.mention}\n"
                    f"**Expires:** {duration}\n"
                    f"**Reason:** {reason}"
                ),
                color=discord.Color.blue(),
            )
            embed.set_thumbnail(url=member_obj.display_avatar.url)
            try:
                return await state.end(success=embed)
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        except Exception as e:
            logger.warning(f"{str(e).capitalize()}")
            raise

    # DONE
    async def handle_vegan_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        duration = None
        is_channel_scope = False
        is_modification = False
        reason = "No reason provided."

        vegans = await Vegan.select(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )
        if not vegans[0]:
            return await state.end(
                warning=f"\U000026a0\U0000fe0f "
                f"{member_obj.mention} is already vegan."
            )

        vegan = Vegan(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )
        await vegan.create()

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=duration,
            is_channel_scope=is_channel_scope,
            is_modification=is_modification,
            member=member_obj,
            message=message,
            reason=reason,
        )
        embed = discord.Embed(
            title=f"\U0001f525\U0001f525 {member_obj.display_name} "
            f"is going Vegan!!!\U0001f525\U0001f525",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Celebrate!** Stick around and do some activism with us!"
            ),
            color=discord.Color.green(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_flag_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        duration = None
        is_channel_scope = False
        is_modification = False
        reason = "No reason provided."
        if len(args) > 1:
            reason = " ".join(args[1:])

        if is_reason_modification and existing_guestroom_alias_event:
            is_modification = True
            modified_reason = " ".join(args[1:]) if len(args) > 1 else ""
            match is_reason_modification:
                case "+":
                    reason = existing_guestroom_alias_event.reason + modified_reason
                case "=" | "-":
                    reason = modified_reason
            set_kwargs = {"reason": reason}
            where_kwargs = {
                "channel_snowflake": channel_obj.id,
                "guild_snowflake": message.guild.id,
                "member_snowflake": member_obj.id,
            }
            await Flag.update(set_kwargs=set_kwargs, where_kwargs=where_kwargs)
        else:
            if not is_modification:
                flags = await Flag.select(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=message.guild.id,
                    member_snowflake=member_obj.id,
                )
                if flags:
                    try:
                        return await state.end(
                            warning=f"\U000026a0\U0000fe0f "
                            f"An existing flag already exists for "
                            f"{member_obj.mention}. Modify the flag by "
                            f"putting a single `+`, `-` or `=` "
                            f"and a reason to update the reason."
                        )
                    except Exception as e:
                        return await state.end(
                            error=f"\u274c " f"{str(e).capitalize()}"
                        )

        flag = Flag(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
            reason=reason,
        )
        await flag.create()

        bot = DiscordBot.get_instance()
        cog = bot.get_cog("EventListeners")
        cog.flags.append(flag)

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=duration,
            is_channel_scope=is_channel_scope,
            is_modification=is_modification,
            member=member_obj,
            message=message,
            reason=reason,
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member_obj.display_name} Flagged",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Reason:** {reason}"
            ),
            color=discord.Color.red(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_role_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        duration = None
        is_channel_scope = False
        is_modification = False
        reason = "No reason provided."

        try:
            role = await resolve_role(
                ctx_interaction_or_message=message, role_str=alias.role_snowflake
            )
        except Exception as e:
            try:
                logger.warning(f"{str(e).capitalize()}")
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f "
                    f"Role `{alias.role_snowflake}` was not found."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            await member_obj.add_roles(role, reason="Added role")
        except discord.Forbidden as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=duration,
            is_channel_scope=is_channel_scope,
            is_modification=is_modification,
            member=member_obj,
            message=message,
            reason=reason,
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member_obj.display_name} Roled",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.blurple(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_text_mute_alias(
        self,
        alias,
        action_information,
        channel,
        member,
        message,
        state
    ):
        if not action_information['action_modification']:
            text_mute = action_information['alias_class'](
                channel_snowflake=action_information['action_channel_snowflake'],
                expires_in=action_information['action_expires_in'],
                guild_snowflake=action_information['action_guild_snowflake'],
                member_snowflake=action_information['action_member_snowflake'],
                reason=action_information['action_reason'],
            )
            await text_mute.create()

        try:
            await channel.set_permissions(
                member=member,
                send_messages=False,
                add_reactions=False,
                reason=action_information['action_reason'],
            )
        except discord.Forbidden as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

        await History.send_entry(
            alias=alias,
            channel=channel,
            duration=action_information['action_duration'],
            is_channel_scope=False,
            is_modification=action_information['action_modification'],
            member=member,
            message=message,
            reason=action_information['action_reason'],
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} Text Muted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {action_information['action_duration']}\n"
                f"**Reason:** {action_information['action_reason']}"
            ),
            color=discord.Color.green(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_voice_mute_alias(
        self,
        alias,
        action_information,
        channel,
        member,
        message,
        state
    ):
        if not action_information['action_modification']:
            voice_mute = action_information['alias_class'](
                channel_snowflake=action_information['action_channel_snowflake'],
                expires_in=action_information['action_expires_in'],
                guild_snowflake=action_information['action_guild_snowflake'],
                member_snowflake=action_information['action_member_snowflake'],
                reason=action_information['action_reason'],
                target='user',
            )
            await voice_mute.create()

        is_channel_scope = False
        if member.voice and member.voice.channel:
            if member.voice.channel.id == action_information['action_channel_snowflake']:
                is_channel_scope = True
                try:
                    await member.edit(mute=True, reason=action_information['action_reason'])
                except discord.Forbidden as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")

        await History.send_entry(
            alias=alias,
            channel=channel,
            duration=action_information['action_duration'],
            is_channel_scope=is_channel_scope,
            is_modification=action_information['action_modification'],
            member=member,
            message=message,
            reason=action_information['action_reason'],
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member.display_name} has been Voice Muted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {action_information['action_duration']}\n"
                f"**Reason:** {action_information['action_reason']}"
            ),
            color=discord.Color.green(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_unban_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        is_channel_scope = False

        bans = await Ban.select(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )
        if not bans[0]:
            return await state.end(
                warning=f"\U000026a0\U0000fe0f "
                f"{member_obj.mention} is not currently banned "
                f"in {channel_obj.mention}."
            )
        if bans[0].expires_in is None and executor_role == "Moderator":
            try:
                return await state.end(
                    warning="\U000026a0\U0000fe0f "
                    "Only coordinators and above can undo permanent bans."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

        await Ban.delete(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )

        try:
            await channel_obj.set_permissions(member_obj, overwrite=None)
        except discord.Forbidden as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

        if member_obj.voice and member_obj.voice.channel:
            if member_obj.voice.channel.id == channel_obj.id:
                is_channel_scope = True
                try:
                    await member_obj.move_to(None, reason="Unbanned")
                except discord.Forbidden as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=None,
            is_channel_scope=is_channel_scope,
            is_modification=True,
            member=member_obj,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member_obj.display_name} has been Unbanned",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_carnist_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):

        vegans = await Vegan.select(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )
        if not vegans[0]:
            return await state.end(
                warning=f"\U000026a0\U0000fe0f "
                f"{member_obj.mention} is not currently vegan."
            )
        
        await Vegan.delete(
            channel_snowflake=channel_obj.id,
            member_snowflake=member_obj.id,
            guild_snowflake=message.guild.id,
        )

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=None,
            is_channel_scope=False,
            is_modification=True,
            member=member_obj,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"\U0001f44e\U0001f44e "
            f"{member_obj.display_name} is a Carnist \U0001f44e\U0001f44e",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.red(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_unflag_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        flags = await Flag.select(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )
        if not flags:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f "
                    f"{member_obj.mention} is not currently flagged in "
                    f"{channel_obj.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

        bot = DiscordBot.get_instance()
        cog = bot.get_cog("EventListeners")
        for flag in cog.flags:
            if flag.channel_snowflake == channel_obj.id:
                cog.flags.remove(flag)
                break

        await Flag.delete(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=None,
            is_channel_scope=False,
            is_modification=True,
            member=member_obj,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member_obj.display_name} has been Unflagged",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_unmute_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        is_channel_scope = False

        voice_mutes = await VoiceMute.select(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
            target="user",
        )
        if not voice_mutes:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f "
                    f"{member_obj.mention} is not currently voice-muted "
                    f"in {channel_obj.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        if voice_mutes[0].expires_in is None:
            if executor_role == "Moderator":
                try:
                    return await state.end(
                        warning="\U000026a0\U0000fe0f "
                        "Only coordinators and above can undo "
                        "permanent voice-mutes."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")

        await VoiceMute.delete(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
            target="user",
        )

        if member_obj.voice and member_obj.voice.channel:
            try:
                is_channel_scope = True
                await member_obj.edit(mute=False)
            except discord.Forbidden as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=None,
            is_channel_scope=is_channel_scope,
            is_modification=True,
            member=member_obj,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member_obj.display_name} has been Unmuted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_unrole_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        try:
            role = await resolve_role(
                ctx_interaction_or_message=message, role_str=alias.role_snowflake
            )
        except Exception as e:
            try:
                logger.warning(f"{str(e).capitalize()}")
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f "
                    f"Role `{alias.role_snowflake}` was not found."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        if role not in member_obj.roles:
            try:
                return await state.end(
                    warning=f"{get_random_emoji()} "
                    f"{member_obj.mention} does not have {role.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        try:
            await member_obj.remove_roles(role)
        except discord.Forbidden as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=None,
            is_channel_scope=False,
            is_modification=True,
            member=member_obj,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member_obj.display_name} has been Unroled",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

    # DONE
    async def handle_untextmute_alias(
        self,
        alias,
        args,
        channel_obj,
        executor_role,
        existing_guestroom_alias_event,
        is_duration_modification,
        is_reason_modification,
        member_obj,
        message,
        state,
    ):
        text_mute = await TextMute.select(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )
        if not text_mute:
            try:
                return await state.end(
                    warning=f"\U000026a0\U0000fe0f "
                    f"{member_obj.mention} is not currently text-muted "
                    f"in {channel_obj.mention}."
                )
            except Exception as e:
                return await state.end(error=f"\u274c {str(e).capitalize()}")
        if text_mute.expires_in is None:
            if executor_role == "Moderator":
                try:
                    return await state.end(
                        warning="\U000026a0\U0000fe0f "
                        "Only coordinators and above can undo "
                        "permanent text-mutes."
                    )
                except Exception as e:
                    return await state.end(error=f"\u274c {str(e).capitalize()}")

        await TextMute.delete(
            channel_snowflake=channel_obj.id,
            guild_snowflake=message.guild.id,
            member_snowflake=member_obj.id,
        )

        try:
            await channel_obj.set_permissions(member=member_obj, send_messages=None)
        except discord.Forbidden as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")

        await History.send_entry(
            alias=alias,
            channel=channel_obj,
            duration=None,
            is_channel_scope=False,
            is_modification=True,
            member=member_obj,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} "
            f"{member_obj.display_name} has been Unmuted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member_obj.mention}\n"
                f"**Channel:** {channel_obj.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member_obj.display_avatar.url)
        try:
            return await state.end(success=embed)
        except Exception as e:
            return await state.end(error=f"\u274c {str(e).capitalize()}")


async def setup(bot: DiscordBot):
    cog = Aliases(bot)
    await bot.add_cog(cog)