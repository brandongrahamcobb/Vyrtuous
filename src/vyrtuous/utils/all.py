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
from datetime import datetime, timezone
from typing import Optional
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.duration import DurationObject
from vyrtuous.utils.emojis import Emojis
from vyrtuous.utils.setup_logging import logger
import discord

class All:
        
    def __init__(self):
        self.bot = DiscordBot.get_instance()

    @classmethod
    async def create_pages_to_show_channels_by_guild_and_member(cls, channel_snowflakes: list[int | None], guild_snowflake: int, member_snowflake: int, member_type):
        if channel_snowflakes is None:
            return None
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        channel_mentions = []
        for channel_snowflake in channel_snowflakes:
            if not channel_snowflake:
                continue
            channel = bot.get_channel(channel_snowflake)
            if not channel:
                logger.warning("Channel not found in create_pages_to_show_channels_by_guild_and_member")
                continue
            channel_mentions.append(channel.mention)
        user = bot.get_user(member_snowflake)
        if not user:
            logger.warning("Member not found in create_pages_to_show_channels_by_guild_and_member")
            return
        pages = []
        chunk_size = 18
        for i in range(0, len(channel_mentions), chunk_size):
            chunk = channel_mentions[i:i+chunk_size]
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {member_type.SINGULAR} Channels for {user.display_name}',
                description = '\n'.join(f'• {channel}' for channel in chunk),
                color = discord.Color.gold()
            )
            pages.append(embed)
        return pages

    @classmethod
    async def create_pages_to_show_guilds_by_member(cls, guild_snowflakes: list[int], member_snowflake: int, member_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        pages = []
        user = bot.get_user(member_snowflake)
        if not user:
            logger.warning("Member not found in create_pages_to_show_channels_by_guild_and_member")
            return
        embed = discord.Embed(
            title = f'{emoji.get_random_emoji()} {user.display_name} as {member_type.SINGULAR} in All Servers',
            description = ', '.join(str(guild_snowflakes)),
            color = discord.Color.blurple()
        )
        pages.append(embed)
        return pages

    @classmethod
    async def create_pages_to_show_members_by_channel_and_guild(cls, channel_snowflake: int, guild_snowflake: int, members, member_type):
        if members is None:
            return None
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        lines = []
        for member in members:
            user = bot.get_user(member.member_snowflake)
            if not user:
                logger.warning("Member not found in create_pages_to_show_members_by_channel_and_guild")
                continue
            lines.append(f'• {member.member_mention}')
        channel = bot.get_channel(channel_snowflake)
        pages = []
        chunk_size = 18
        for i in range(0, len(lines), chunk_size):
            chunk = lines[i:i+chunk_size]
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {member_type.PLURAL} in {channel.mention}',
                description='\n'.join(chunk),
                color=discord.Color.gold()
            )
            pages.append(embed)
        return pages
    

    @classmethod
    async def create_pages_to_show_members_by_guild(cls, guild_snowflake: int, members, member_type):
        if members is None:
            return None
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        lines = []
        for member in members:
            user = bot.get_user(member.member_snowflake)
            if not user:
                logger.warning("Member not found in create_pages_to_show_members_by_guild")
                continue
            lines.append(f'• {member.member_mention}')
        guild = bot.get_guild(guild_snowflake)
        pages = []
        chunk_size = 18
        for i in range(0, len(lines), chunk_size):
            chunk = lines[i:i+chunk_size]
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {member_type.PLURAL} in {guild.name}',
                description='\n'.join(chunk),
                color=discord.Color.gold()
            )
            pages.append(embed)
        return pages
    

    @classmethod
    async def create_pages_to_show_channels_by_guild_and_members(cls, channel_snowflakes: list[int | None], guild_snowflake: int, members, member_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        pages = []
        chunk_size = 18
        guild = bot.get_channel(guild_snowflake)
        for member in members:
            user = bot.get_user(member.member_snowflake)
            if not user:
                logger.warning("Member not found in create_pages_to_show_channels_by_guild_and_members")
                continue
            channel_mentions = []
            for channel_snowflake in channel_snowflakes:
                if not channel_snowflake:
                    logger.warning("Channel not found in create_pages_to_show_channels_by_guild_and_members")
                    continue
                channel = bot.get_channel(channel_snowflake)
                channel_mentions.append(channel.mention)
            for i in range(0, len(channel_mentions), chunk_size):
                chunk = channel_mentions[i:i + chunk_size]
                embed = discord.Embed(
                    title=f'{emoji.get_random_emoji()} {member_type.SINGULAR} Channels for {user.display_name} in {guild.name}',
                    description='\n'.join(f'• {channel}' for channel in chunk),
                    color=discord.Color.gold()
                )
                pages.append(embed)
        return pages

    @classmethod
    async def create_pages_to_show_guilds_by_members(cls, members, member_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        pages = []
        for member in members:
            user = bot.get_user(member.member_snowflake)
            if not user:
                logger.warning("Member not found in create_pages_to_show_guilds_by_members")
                continue
            guild = bot.get_guild(member.guild_snowflake)
            if not guild:
                continue
            embed = discord.Embed(
                title=f'{emoji.get_random_emoji()} {user.display_name} is a {member_type.SINGULAR} in {guild.name}',
                description=f'• {guild.name}',
                color=discord.Color.blurple()
            )
            pages.append(embed)
        return pages
    
    @classmethod
    async def create_pages_from_moderations_by_channel_and_guild(cls, guild_snowflake, moderations, moderation_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        guild = bot.get_guild(guild_snowflake)
        lines_by_channel = {}
        pages = []
        for moderation in moderations:
            channel = bot.get_channel(moderation.channel_snowflake)
            user = bot.get_user(moderation.member_snowflake)
            if not channel or not user:
                continue
            if moderation.expires_at is None:
                duration = DurationObject(0)
            else:
                delta = moderation.expires_at - datetime.now(timezone.utc)
                duration = DurationObject.from_timedelta(delta)
            lines_by_channel.setdefault(channel.mention, []).append(f'{user.mention}\nReason: {moderation.reason}\nDuration: {str(duration)}')
        for channel_mention, entries in lines_by_channel.items():
            chunk_size = 18
            for i in range(0, len(entries), chunk_size):
                embed = discord.Embed(
                    title=f'{emoji.get_random_emoji()} {moderation_type.PLURAL} for {channel_mention} in {guild.name}',
                    color=discord.Color.red()
                )
                embed.add_field(name='Users', value='\n\n'.join(entries[i:i + chunk_size]), inline=False)
                pages.append(embed)
        return pages

    @classmethod
    async def create_pages_from_moderations_by_guild_and_member(cls, guild_snowflake, moderations, moderation_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        guild = bot.get_guild(guild_snowflake)
        lines_by_user = {}
        pages = []
        for moderation in moderations:
            channel = bot.get_channel(moderation.channel_snowflake)
            user = bot.get_user(moderation.member_snowflake)
            if not channel or not user:
                continue
            if moderation.expires_at is None:
                duration = DurationObject(0)
            else: 
                delta = moderation.expires_at - datetime.now(timezone.utc)
                duration = DurationObject.from_timedelta(delta)
            lines_by_user.setdefault(user, []).append(
                f'{channel.mention}\nReason: {moderation.reason}\nDuration: {str(duration)}'
            )
        for user, entries in lines_by_user.items():
            chunk_size = 18
            for i in range(0, len(entries), chunk_size):
                embed = discord.Embed(
                    title=f'{emoji.get_random_emoji()} {moderation_type.PLURAL} for {user.display_name} in {guild.name}',
                    color=discord.Color.red()
                )
                embed.add_field(name='Channels', value='\n\n'.join(entries[i:i + chunk_size]), inline=False)
                pages.append(embed)
        return pages


    @classmethod
    async def create_pages_from_moderations_by_guild(cls, guild_snowflake, moderations, moderation_type):
        bot = DiscordBot.get_instance()
        emoji = Emojis()
        guild = bot.get_guild(guild_snowflake)
        lines_by_guild = []
        pages = []
        for moderation in moderations:
            user = bot.get_user(moderation.member_snowflake)
            if not user:
                continue
            if moderation.expires_at is None:
                duration = DurationObject(0)
            else: 
                delta = moderation.expires_at - datetime.now(timezone.utc)
                duration = DurationObject.from_timedelta(delta)
            lines_by_guild.append(f'• {user.mention} — {str(duration)}')
        chunk_size = 18
        for i in range(0, len(lines_by_guild), chunk_size):
            chunk = lines_by_guild[i:i+chunk_size]
            embed = discord.Embed(
                title = f'{emoji.get_random_emoji()} {moderation_type.PLURAL} in {guild.name}',
                description = '\n'.join(chunk),
                color = discord.Color.red()
            )
            pages.append(embed)
        return pages