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
    async def create_pages_to_show_channels(cls, channel_snowflakes: list[int | None], member_type):
        if channel_snowflakes is None:
            return None
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        channel_mentions = []
        for channel_snowflake in channel_snowflakes:
            if not channel_snowflake:
                continue
            channel = bot.get_channel(channel_snowflake)
            channel_mentions.append(channel.mention if vc else f'Unknown Channel ({channel_snowflake})')
        pages = []
        chunk_size = 18
        for i in range(0, len(channel_mentions), chunk_size):
            chunk = channel_mentions[i:i+chunk_size]
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {member_type.SINGULAR}',
                description = '\n'.join(f'• {channel}' for channel in chunk),
                color = discord.Color.gold()
            )
            pages.append(embed)
        return pages

    @classmethod
    async def create_pages_to_show_guilds_by_members(cls, members, member_type):
        if members is None:
            return None
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        pages = []
        for member in members:
            user = bot.get_user(member.member_snowflake)
            name = user.name if user else f'User ID {member.member_snowflake}'
            embed = discord.Embed(
                title = f'{emoji.get_random_emoji()} {member_type.PLURAL}',
                description = ', '.join(str(member.guild_snowflakes)) if member.guild_snowflakes else 'No known guilds',
                color = discord.Color.blurple()
            )
            pages.append(embed)
        return pages

    @classmethod
    async def create_pages_to_show_members(cls, members, member_type):
        if members is None:
            return None
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        lines = []
        for member in members:
            user = bot.get_user(member.member_snowflake)
            if user:
                lines.append(f'• {user.display_name} — <@{member.member_snowflake}>')
        pages = []
        chunk_size = 18
        for i in range(0, len(lines), chunk_size):
            chunk = lines[i:i+chunk_size]
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {member_type.PLURAL}',
                description='\n'.join(chunk),
                color=discord.Color.gold()
            )
            pages.append(embed)
        return pages

    @classmethod
    async def create_pages_to_show_guilds_by_member(cls, guilds, member_snowflake, member_type):
        if guilds is None:
            return None
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        pages = []
        user = bot.get_user(member_snowflake)
        name = user.name if user else f'User ID {member_snowflake}'
        embed = discord.Embed(
            title = f'{emoji.get_random_emoji()} {member_type.PLURAL}',
            description = ', '.join(guilds) if guilds else 'No known guilds',
            color = discord.Color.blurple()
        )
        pages.append(embed)
        return pages