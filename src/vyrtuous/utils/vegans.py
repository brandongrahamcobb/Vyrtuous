''' vegans.py A utility module for managing vegan members in the Vyrtuous Discord bot.
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

import discord

class Vegans:

    state: bool = False
    vegans = set()

    @classmethod
    async def unrestrict(cls, guild: discord.Guild, member: discord.Member):
        bot = DiscordBot.get_instance()
        uid = member.id
        async with bot.db_pool.acquire() as conn:
            ban_rows = await conn.fetch('SELECT channel_id FROM active_bans WHERE discord_snowflake=$1', uid)
            mute_rows = await conn.fetch('SELECT channel_id FROM active_voice_mutes WHERE discord_snowflake=$1', uid)
            text_rows = await conn.fetch('SELECT channel_id FROM active_text_mutes WHERE discord_snowflake=$1', uid)
        for r in ban_rows:
            try: await guild.unban(discord.Object(id=uid), reason='Toggle OFF')
            except: pass
        for r in mute_rows:
            ch = guild.get_channel(r['channel_id'])
            if ch and member.voice and member.voice.mute: await member.edit(mute=False)
        for r in text_rows:
            ch = guild.get_channel(r['channel_id'])
            text_mute_role = discord.utils.get(guild.roles, name='TextMuted')
            if ch and text_mute_role and text_mute_role in member.roles: await member.remove_roles(text_mute_role)
        async with bot.db_pool.acquire() as conn:
            await conn.execute('DELETE FROM active_bans WHERE discord_snowflake=$1', uid)
            await conn.execute('DELETE FROM active_voice_mutes WHERE discord_snowflake=$1', uid)
            await conn.execute('DELETE FROM active_text_mutes WHERE discord_snowflake=$1', uid)

    @classmethod
    def add_vegan(cls, member_id: int):
        cls.vegans.add(member_id)
        
    @classmethod
    def get_vegans(cls):
        return cls.vegans

    @classmethod
    def remove_vegan(cls, member_id: int):
        cls.vegans.discard(member_id)
        
    @classmethod
    def toggle_state(cls):
        cls.state = not cls.state
        return cls.state


