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

from vyrtuous.db.actions.ban import Ban
from vyrtuous.db.actions.flag import Flag
from vyrtuous.db.actions.hide import Hide
from vyrtuous.db.actions.role import Role
from vyrtuous.db.actions.text_mute import TextMute
from vyrtuous.db.roles.vegan import Vegan
from vyrtuous.db.actions.voice_mute import VoiceMute
from vyrtuous.db.mgmt.stream import Streaming
from vyrtuous.utils.invincibility import Invincibility


class Aliases(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.alias_help = {
            "ban": [
                "**member** (Required): Tag a member or include their ID",
                "**duration** (Optional): #(m|h|d)\n0 = permanent",
                "**reason** (Optional): Reason for ban",
            ],
            "flag": [
                "**member** (Required): Tag a member or include their ID",
                "**reason** (Optional): Reason for flag",
            ],
            "hide": [
                "**member** (Required): Tag a member or include their ID",
                "**duration** (Optional): #(m|h|d)\n0 = permanent",
                "**reason** (Optional): Reason for hiding",
            ],
            "role": [
                "**member** (Required): Tag a member or include their ID",
                "**role** (Required): Tag a role or include its ID",
            ],
            "tmute": [
                "**member** (Required): Tag a member or include their ID",
                "**duration** (Optional): #(m|h|d)\n0 = permanent",
                "**reason** (Required): Reason for text-mute",
            ],
            "vegan": ["**member** (Required): Tag a member or include their ID"],
            "vmute": [
                "**member** (Required): Tag a member or include their ID",
                "**duration** (Optional): #(m|h|d)\n0 = permanent",
                "**reason** (Optional): Reason for voice-mute",
            ],
        }
        self.category_to_description = {
            "ban": "Toggles a ban.",
            "flag": "Toggles a moderation flag.",
            "hide": "Toggles a hidden channel.",
            "role": "Toggles a role to a user.",
            "tmute": "Toggles a mute in text channels.",
            "vegan": "Toggles a vegan.",
            "vmute": "Toggles a mute in voice channels.",
        }
        self.category_to_permission_level = {
            "ban": "Moderator",
            "flag": "Moderator",
            "hide": "Coordinator",
            "role": "Coordinator",
            "tmute": "Moderator",
            "vegan": "Moderator",
            "vmute": "Moderator",
        }
        self.bot = bot
        self.invincible_members = Invincibility.get_invincible_members()

    async def handle_ban_alias(self, alias, action_information, member, message, state):

        ban = Ban(
            channel_snowflake=action_information["action_channel_snowflake"],
            expires_in=action_information["action_expires_in"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
            reason=action_information["action_reason"],
        )
        await ban.create()

        role = message.guild.get_role(action_information["action_role_snowflake"])
        if not role:
            return await state.end(
                warning=f"Role `{action_information['action_role_snowflake']}` was not found."
            )

        await Role.administer_role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )

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
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=action_information["action_duration"],
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
        )

        embed = await Ban.act_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    async def handle_unban_alias(
        self, alias, action_information, member, message, state
    ):

        await Ban.delete(
            channel_snowflake=action_information["action_channel_snowflake"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
        )
        is_channel_scope = False

        await Role.revoke_role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )

        if member.voice and member.voice.channel:
            try:
                is_channel_scope = True
                await member.edit(mute=False)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=None,
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
        )

        embed = await Ban.undo_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    async def handle_hide_alias(
        self, alias, action_information, member, message, state
    ):

        hide = Hide(
            channel_snowflake=action_information["action_channel_snowflake"],
            expires_in=action_information["action_expires_in"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
            reason=action_information["action_reason"],
        )
        await hide.create()

        role = message.guild.get_role(action_information["action_role_snowflake"])
        if not role:
            return await state.end(
                warning=f"Role `{action_information['action_role_snowflake']}` was not found."
            )

        await Role.administer_role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )

        is_channel_scope = False
        if (
            member.voice
            and member.voice.channel
            and member.voice.channel.id
            == action_information["action_channel_snowflake"]
        ):
            is_channel_scope = True
            try:
                await member.move_to(None, reason=action_information["action_reason"])
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=action_information["action_duration"],
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
        )

        embed = await Hide.act_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_vegan_alias(
        self, alias, action_information, member, message, state
    ):

        vegan = Vegan(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
        )
        await vegan.create()

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provied.",
        )

        embed = await Vegan.act_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_flag_alias(
        self, alias, action_information, member, message, state
    ):

        flag = Flag(
            channel_snowflake=action_information["action_channel_snowflake"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            reason=action_information["action_reason"],
        )
        await flag.create()

        bot = DiscordBot.get_instance()
        cog = bot.get_cog("ChannelEventListeners")
        cog.flags.append(flag)

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=action_information["action_duration"],
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
        )

        embed = await Flag.act_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_role_alias(
        self, alias, action_information, member, message, state
    ):
        db_role = Role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )
        await db_role.create()

        role = message.guild.get_role(alias.role_snowflake)
        if not role:
            return await state.end(
                warning=f"Role `{alias.role_snowflake}` was not found."
            )

        await Role.administer_role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=action_information["action_duration"],
            is_channel_scope=False,
            is_modification=False,
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await Role.act_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_text_mute_alias(
        self, alias, action_information, member, message, state
    ):

        text_mute = TextMute(
            channel_snowflake=action_information["action_channel_snowflake"],
            expires_in=action_information["action_expires_in"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
            reason=action_information["action_reason"],
        )
        await text_mute.create()

        await Role.administer_role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=action_information["action_duration"],
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
        )

        embed = await TextMute.act_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_voice_mute_alias(
        self, alias, action_information, member, message, state
    ):

        voice_mute = VoiceMute(
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
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=action_information["action_duration"],
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason=action_information["action_reason"],
        )

        embed = await VoiceMute.act_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_unhide_alias(
        self, alias, action_information, member, message, state
    ):
        await Hide.delete(
            channel_snowflake=action_information["action_channel_snowflake"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
        )
        is_channel_scope = False

        await Role.revoke_role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=None,
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await Hide.undo_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_carnist_alias(
        self, alias, action_information, member, message, state
    ):
        await Vegan.delete(
            channel_snowflake=action_information["action_channel_snowflake"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
        )
        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await Vegan.undo_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_unflag_alias(
        self, alias, action_information, member, message, state
    ):
        await Flag.delete(
            channel_snowflake=action_information["action_channel_snowflake"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
        )
        bot = DiscordBot.get_instance()
        cog = bot.get_cog("ChannelEventListeners")
        for flag in cog.flags:
            if flag.channel_snowflake == action_information["action_channel_snowflake"]:
                cog.flags.remove(flag)
                break

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await Flag.undo_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_unmute_alias(
        self, alias, action_information, member, message, state
    ):
        await VoiceMute.delete(
            channel_snowflake=action_information["action_channel_snowflake"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
        )
        is_channel_scope = False

        if member.voice and member.voice.channel:
            try:
                is_channel_scope = True
                await member.edit(mute=False)
            except discord.Forbidden as e:
                return await state.end(error=str(e).capitalize())

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=None,
            is_channel_scope=is_channel_scope,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await VoiceMute.undo_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_unrole_alias(
        self, alias, action_information, member, message, state
    ):
        await Role.delete(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )

        role = message.guild.get_role(alias.role_snowflake)
        if not role:
            return await state.end(
                warning=f"Role `{alias.role_snowflake}` was not found."
            )

        await Role.revoke_role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
        )

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await Role.undo_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)

    # DONE
    async def handle_untextmute_alias(
        self, alias, action_information, member, message, state
    ):
        await TextMute.delete(
            channel_snowflake=action_information["action_channel_snowflake"],
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
        )

        await Role.revoke_role(
            guild_snowflake=action_information["action_guild_snowflake"],
            member_snowflake=action_information["action_member_snowflake"],
            role_snowflake=action_information["action_role_snowflake"],
            state=state,
        )

        await Streaming.send_entry(
            alias=alias,
            channel_snowflake=action_information["action_channel_snowflake"],
            duration=None,
            is_channel_scope=False,
            is_modification=action_information["action_modification"],
            member=member,
            message=message,
            reason="No reason provided.",
        )

        embed = await TextMute.undo_embed(
            action_information=action_information, source=message
        )

        return await state.end(success=embed)


async def setup(bot: DiscordBot):
    cog = Aliases(bot)
    await bot.add_cog(cog)
