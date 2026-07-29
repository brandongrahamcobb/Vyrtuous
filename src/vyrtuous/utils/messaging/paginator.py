"""!/bin/python3
paginator.py  The purpose of this program is to provide a pagintor for scrolling through embeds.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import asyncio
from typing import Union

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.errors.error import ChannelNotFound, CheckFailure, GuildNotFound


class Paginator:
    NAV_EMOJIS = {"\u2b05\ufe0f": -1, "\u27a1\ufe0f": 1}

    def __init__(self):
        self.current_page = 0
        self._reaction_lock = asyncio.Lock()

    async def start_by_message(
        self,
        pages: list[discord.Embed],
        source: Union[commands.Context, discord.Interaction, discord.Message],
        *,
        ephemeral=False,
        timeout=60,
    ) -> discord.Message:
        bot: DiscordBot = DiscordBot.get_instance()
        if source.guild is None:
            raise CheckFailure("This message must be sent in a guild.")
        embed = self.get_current_embed(guild_snowflake=source.guild.id, pages=pages)
        if isinstance(source, discord.Interaction):
            if not source.response.is_done():
                await source.response.send_message(embed=embed, ephemeral=ephemeral)
            message: discord.Message | discord.interactions.InteractionMessage = (
                await source.original_response()
            )
        else:
            message = await source.reply(embed=embed)
        for emoji in self.NAV_EMOJIS:
            await message.add_reaction(emoji)
        bot.loop.create_task(
            self.wait_for_reactions(message=message, pages=pages, timeout=timeout)
        )
        return message

    async def start_without_message(
        self,
        channel_snowflake: int,
        guild_snowflake: int,
        pages: list[discord.Embed],
        *,
        timeout=60,
    ) -> discord.Message:
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise GuildNotFound(str(guild_snowflake))
        channel = guild.get_channel(channel_snowflake)
        if channel is None:
            raise ChannelNotFound(str(channel_snowflake))
        if not isinstance(
            channel,
            (discord.TextChannel, discord.VoiceChannel, discord.StageChannel),
        ):
            raise commands.CheckFailure("This message must be sent to a valid channel.")
        embed = self.get_current_embed(guild_snowflake=guild_snowflake, pages=pages)
        message = await channel.send(embed=embed)
        for emoji in self.NAV_EMOJIS:
            await message.add_reaction(emoji)
        bot.loop.create_task(
            self.wait_for_reactions(message=message, pages=pages, timeout=timeout)
        )
        return message

    def get_current_embed(
        self, guild_snowflake: int, pages: list[discord.Embed]
    ) -> discord.Embed:
        bot: DiscordBot = DiscordBot.get_instance()
        guild = bot.get_guild(guild_snowflake)
        if guild is None:
            raise GuildNotFound(str(guild_snowflake))
        embed = pages[self.current_page].copy()
        total_pages = len(pages)
        label = "page"
        embed.set_footer(
            text=f"{label} {self.current_page + 1}/{total_pages} • {guild.name}"
        )
        return embed

    async def wait_for_reactions(
        self,
        message: discord.Message,
        pages: list[discord.Embed],
        timeout: int,
    ) -> None:
        bot: DiscordBot = DiscordBot.get_instance()

        def look(reaction, user):
            return (
                reaction.message.id == message.id
                and str(reaction.emoji) in self.NAV_EMOJIS
                and not user.bot
            )

        while True:
            try:
                reaction, user = await bot.wait_for(
                    "reaction_add", timeout=timeout, check=look
                )
            except asyncio.TimeoutError:
                try:
                    await message.clear_reactions()
                except Exception as e:
                    bot.logger.warning(str(e).capitalize())
                break
            await self.handle_reaction(
                message=message,
                pages=pages,
                reaction=reaction,
            )
            try:
                await message.remove_reaction(reaction.emoji, user)
            except Exception as e:
                bot.logger.warning(str(e).capitalize())

    async def handle_reaction(self, message, pages, reaction) -> None:
        async with self._reaction_lock:
            action = self.NAV_EMOJIS[str(reaction.emoji)]
            if isinstance(action, int):
                self.current_page = max(
                    0, min(self.current_page + action, len(pages) - 1)
                )
                await message.edit(
                    embed=self.get_current_embed(
                        guild_snowflake=message.channel.guild.id,
                        pages=pages,
                    )
                )
