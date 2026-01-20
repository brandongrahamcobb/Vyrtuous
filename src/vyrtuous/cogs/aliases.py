"""aliases.py A discord.py cog containing command aliases for the Vyrtuous bot.

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

from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.database.settings.streaming import Streaming

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
        self, alias, action_information, channel, member, message, state
    ):
        if not action_information["action_modification"]:
            ban = action_information["alias_class"](
                channel_snowflake=action_information["action_channel_snowflake"],
                expires_in=action_information["action_expires_in"],
                guild_snowflake=action_information["action_guild_snowflake"],
                member_snowflake=action_information["action_member_snowflake"],
                reason=action_information["action_reason"],
            )
            await ban.create()

        try:
            await channel.set_permissions(
                member, view_channel=False, reason=action_information["action_reason"]
            )
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())
        if (
            member.voice
            and member.voice.channel
            and member.voice.channel.id == channel.id
        ):
            is_channel_scope = True
            try:
                await member.move_to(None, reason=action_information["action_reason"])
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=action_information["action_duration"],
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Banned",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Expires:** {action_information["action_duration"]}\n"
                f"**Reason:** {action_information["action_modification"]}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_vegan_alias(
        self, alias, action_information, channel, member, message, state
    ):
        if not action_information["action_modification"]:
            vegan = action_information["alias_class"](
                guild_snowflake=action_information["action_guild_snowflake"],
                member_snowflake=action_information["action_member_snowflake"],
            )
            await vegan.create()

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provied.",
        )
        embed = discord.Embed(
            title=f"\U0001f525\U0001f525 {member.display_name} "
            f"is going Vegan!!!\U0001f525\U0001f525",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Celebrate!** Stick around and do some activism with us!"
            ),
            color=discord.Color.green(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_flag_alias(
        self, alias, action_information, channel, member, message, state
    ):
        if not action_information["action_modification"]:
            flag = action_information["alias_class"](
                channel_snowflake=action_information["action_channel_snowflake"],
                guild_snowflake=action_information["action_guild_snowflake"],
                member_snowflake=action_information["action_member_snowflake"],
                reason=action_information["action_reason"],
            )
            await flag.create()

        bot = DiscordBot.get_instance()
        cog = bot.get_cog("EventListeners")
        cog.flags.append(flag)

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=action_information["action_duration"],
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} Flagged",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Reason:** {action_information['action_reason']}"
            ),
            color=discord.Color.red(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_role_alias(
        self, alias, action_information, channel, member, message, state
    ):
        if not action_information["action_modification"]:
            role_obj = action_information["alias_class"](
                guild_snowflake=action_information["action_guild_snowflake"],
                member_snowflake=action_information["action_member_snowflake"],
                role_snowflake=action_information["action_role_snowflake"],
            )
            await role_obj.create()
        role = message.guild.get_role(alias.role_snowflake)
        if not role:
            return await state.end(
                warning=f"Role `{alias.role_snowflake}` was not found."
            )
        try:
            await member.add_roles(role, reason="Added role")
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=action_information["action_duration"],
            is_channel_scope=False,
            is_modification=False,
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} Roled",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.blurple(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_text_mute_alias(
        self, alias, action_information, channel, member, message, state
    ):
        if not action_information["action_modification"]:
            text_mute = action_information["alias_class"](
                channel_snowflake=action_information["action_channel_snowflake"],
                expires_in=action_information["action_expires_in"],
                guild_snowflake=action_information["action_guild_snowflake"],
                member_snowflake=action_information["action_member_snowflake"],
                reason=action_information["action_reason"],
            )
            await text_mute.create()

        try:
            await channel.set_permissions(
                target=member,
                send_messages=False,
                add_reactions=False,
                reason=action_information["action_reason"],
            )
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=action_information["action_duration"],
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
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

        return await state.end(success=embed)

    # DONE
    async def handle_voice_mute_alias(
        self, alias, action_information, channel, member, message, state
    ):
        if not action_information["action_modification"]:
            voice_mute = action_information["alias_class"](
                channel_snowflake=action_information["action_channel_snowflake"],
                expires_in=action_information["action_expires_in"],
                guild_snowflake=action_information["action_guild_snowflake"],
                member_snowflake=action_information["action_member_snowflake"],
                reason=action_information["action_reason"],
                target="user",
            )
            await voice_mute.create()

        is_channel_scope = False
        if member.voice and member.voice.channel:
            if (
                member.voice.channel.id
                == action_information["action_channel_snowflake"]
            ):
                is_channel_scope = True
                try:
                    await member.edit(
                        mute=True, reason=action_information["action_reason"]
                    )
                except discord.Forbidden as e:
                    return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=action_information["action_duration"],
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
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

        return await state.end(success=embed)

    # DONE
    async def handle_unban_alias(
        self, alias, action_information, channel, member, message, state
    ):
        is_channel_scope = False

        try:
            await channel.set_permissions(member, overwrite=None)
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

        if member.voice and member.voice.channel:
            if member.voice.channel.id == channel.id:
                is_channel_scope = True
                try:
                    await member.move_to(None, reason="Unbanned")
                except discord.Forbidden as e:
                    return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=None,
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Unbanned",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_carnist_alias(
        self, alias, action_information, channel, member, message, state
    ):

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"\U0001f44e\U0001f44e "
            f"{member.display_name} is a Carnist \U0001f44e\U0001f44e",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.red(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_unflag_alias(
        self, alias, action_information, channel, member, message, state
    ):

        bot = DiscordBot.get_instance()
        cog = bot.get_cog("EventListeners")
        for flag in cog.flags:
            if flag.channel_snowflake == channel.id:
                cog.flags.remove(flag)
                break

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Unflagged",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_unmute_alias(
        self, alias, action_information, channel, member, message, state
    ):
        is_channel_scope = False

        if member.voice and member.voice.channel:
            try:
                is_channel_scope = True
                await member.edit(mute=False)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=None,
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Unmuted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_unrole_alias(
        self, alias, action_information, channel, member, message, state
    ):
        role = message.guild.get_role(alias.role_snowflake)
        if not role:
            return await state.end(
                warning=f"Role `{alias.role_snowflake}` was not found."
            )
        if role not in member.roles:
            return await state.end(
                warning=f"{get_random_emoji()} "
                f"{member.mention} does not have {role.mention}."
            )
        try:
            await member.remove_roles(role)
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Unroled",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Role:** {role.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)

    # DONE
    async def handle_untextmute_alias(
        self, alias, action_information, channel, member, message, state
    ):
        try:
            await channel.set_permissions(target=member, send_messages=None)
        except discord.Forbidden as e:
            return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel=channel,
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = discord.Embed(
            title=f"{get_random_emoji()} " f"{member.display_name} has been Unmuted",
            description=(
                f"**By:** {message.author.mention}\n"
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}"
            ),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)

        return await state.end(success=embed)


async def setup(bot: DiscordBot):
    cog = Aliases(bot)
    await bot.add_cog(cog)
