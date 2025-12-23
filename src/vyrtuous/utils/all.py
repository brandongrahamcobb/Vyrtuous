''' temporary_rooms.py A utility module for managing temporary rooms in the Vyrtuous Discord bot.
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
from collections import defaultdict
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.emojis import Emojis
import discord

class All:
        
    def __init__(self):
        self.bot = DiscordBot.get_instance()

    @classmethod
    async def create_show_all_members_pages(cls, guild_name: Optional[str], members, member_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        channel_map = defaultdict(list)
        for member in members:
            channel_map[member.channel_id].append(member.member_id)
        pages = []
        for ch_id, user_ids in sorted(channel_map.items()):
            vc = bot.get_channel(ch_id)
            vc_name = vc.mention if vc else f'Unknown Channel ({ch_id})'
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {member_type.PLURAL} for {vc_name}',
                color=discord.Color.gold())
            for uid in user_ids:
                m = bot.get_user(uid)
                name = m.display_name if m else f'User ID {uid}'
                embed.add_field(name=f'{guild_name}', value=f'• {name} (<@{uid}>)', inline=False)
            pages.append(embed)
        return pages


    @classmethod
    async def create_show_all_channels_pages(cls, member_channel_ids: list[str | None], member_name: Optional[str], member_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        channel_mentions = []
        for ch_id in member_channel_ids:
            if not ch_id:
                continue
            vc = bot.get_channel(ch_id)
            channel_mentions.append(vc.mention if vc else f'Unknown Channel ({ch_id})')
        pages = []
        chunk_size = 18
        for i in range(0, len(channel_mentions), chunk_size):
            chunk = channel_mentions[i:i+chunk_size]
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {member_name} is a {member_type.SINGULAR} in:',
                description = '\n'.join(f'• {ch}' for ch in chunk),
                color = discord.Color.gold()
            )
            pages.append(embed)
        return pages

    @classmethod
    async def create_show_all_members_in_channel_pages(cls, channel_name: Optional[str], member_ids: list[str | None], member_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        lines = []
        for member_id in member_ids:
            m = bot.get_user(member_id)
            if m:
                lines.append(f'• {m.display_name} — <@{member_id}>')
        pages = []
        chunk_size = 18
        for i in range(0, len(lines), chunk_size):
            chunk = lines[i:i+chunk_size]
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {member_type.PLURAL} for {channel_name}',
                description='\n'.join(chunk),
                color=discord.Color.gold()
            )
            pages.append(embed)
        return pages