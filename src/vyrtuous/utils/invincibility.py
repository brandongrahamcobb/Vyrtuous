''' invincibility.py A utility module for granting invincibility to moderation events and revoking all moderation events for a user.

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
'''
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.moderation_action.ban import Ban
from vyrtuous.moderation_action.text_mute import TextMute
from vyrtuous.moderation_action.voice_mute import VoiceMute
import discord

class Invincibility:

    state: bool = False
    invincible_members = set()

    @classmethod
    async def unrestrict(cls, guild_snowflake, member_snowflake):
        target = 'user'
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        member = guild.get_member(member_snowflake)
        bans = await Ban.fetch_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
        text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
        voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake, target='user')
        if bans:
            for ban in bans:
                channel = guild.get_channel(ban.channel_snowflake)
                if channel:
                    try:
                        await channel.set_permissions(member, overwrite=None)
                    except discord.Forbidden:
                        logger.warning(f'Unable to unban {member.name} ({member.id}) in {channel.name} ({channel.id}).')
        if text_mutes:
            for text_mute in text_mutes:
                channel = guild.get_channel(text_mute.channel_snowflake)
                if channel:
                    try:
                        await channel.set_permissions(member, send_messages=True)
                    except discord.Forbidden:
                        logger.warning(f'Unable to untmute {member.name} ({member.id}) in {channel.name} ({channel.id}).')
        if voice_mutes:
            for voice_mute in voice_mutes:
                channel = guild.get_channel(voice_mute.channel_snowflake)
                if channel and member.voice and member.voice.mute:
                    await member.edit(mute=False)
        await Ban.delete_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
        await TextMute.delete_by_guild_and_member(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake)
        await VoiceMute.delete_by_guild_member_and_target(guild_snowflake=guild_snowflake, member_snowflake=member_snowflake, target=target)

    @classmethod
    def add_invincible_member(cls, member_snowflake: int):
        cls.invincible_members.add(member_snowflake)
        
    @classmethod
    def get_invincible_members(cls):
        return cls.invincible_members

    @classmethod
    def remove_invincible_member(cls, member_snowflake: int):
        cls.invincible_members.discard(member_snowflake)
        
    @classmethod
    def toggle_enabled(cls):
        cls.state = not cls.state
        return cls.state


