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
from vyrtuous.utils.ban import Ban
from vyrtuous.utils.text_mute import TextMute
from vyrtuous.utils.voice_mute import VoiceMute
import discord

class Invincibility:

    state: bool = False
    invincible_members = set()

    @classmethod
    async def unrestrict(cls, guild: discord.Guild, member: discord.Member):
        bans = await Ban.fetch_by_guild_and_member(guild_snowflake=guild.id, member_snowflake=member.id)
        text_mutes = await TextMute.fetch_by_guild_and_member(guild_snowflake=guild.id, member_snowflake=member.id)
        voice_mutes = await VoiceMute.fetch_by_guild_member_and_target(guild_snowflake=guild.id, member_snowflake=member.id, target="user")
        if bans:
            for ban in bans:
                try:
                    await guild.unban(ban.member_snowflake, reason='Toggle bans')
                except:
                    pass
        if text_mutes:
            for text_mute in text_mutes:
                ch = guild.get_channel(text_mute.channel_snowflake)
                text_mute_role = discord.utils.get(guild.roles, name='Toggle text mutes')
                if ch and text_mute_role and text_mute_role in member.roles:
                    await member.remove_roles(text_mute_role)
        if voice_mutes:
            for voice_mute in voice_mutes:
                ch = guild.get_channel(voice_mute.channel_snowflake)
                if ch and member.voice and member.voice.mute:
                    await member.edit(mute=False)
        await Ban.delete_by_guild_and_member(guild_snowflake=guild.id, member_snowflake=member.id)
        await TextMute.delete_by_guild_and_member(guild_snowflake=guild.id, member_snowflake=member.id)
        await VoiceMute.delete_by_guild_and_member(guild_snowflake=guild.id, member_snowflake=member.id)

    @classmethod
    def add_invincible_member(cls, member_id: int):
        cls.invincible_members.add(member_id)
        
    @classmethod
    def get_invincible_members(cls):
        return cls.invincible_members

    @classmethod
    def remove_invincible_member(cls, member_id: int):
        cls.invincible_members.discard(member_id)
        
    @classmethod
    def toggle_enabled(cls):
        cls.state = not cls.state
        return cls.state


