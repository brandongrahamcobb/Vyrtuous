''' discord_message_service.py  The purpose of this program is to handle messages in Discord.

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
from discord.ext import commands
from cachetools import TTLCache
from datetime import datetime, timedelta, timezone
from discord.ext import commands
from discord.ui import Button, View
from typing import Union
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.paginator import Paginator
from vyrtuous.utils.developer import Developer
from vyrtuous.utils.setup_logging import logger
from vyrtuous.utils.snowflake import *
from vyrtuous.utils.time_to_complete import TimeToComplete
import discord

cache = TTLCache(maxsize=500, ttl=8*60*60)

class Information(View):
    COLOR_MAP = {
        "green": 0x57F287,
        "yellow": 0xFEE65C,
        "red": 0xED4245
    }

    def __init__(self, message: discord.Message = None):
        super().__init__(timeout=None)
        self.bot = DiscordBot.get_instance()
        self.handler = MessageService(self.bot, self.bot.db_pool)
        self.message_snowflake = message.id
        self._last_called = {}
        self.cooldown_seconds = 10
        self._reported_users = set()

    @discord.ui.button(label="Information", style=discord.ButtonStyle.secondary)
    async def inform(self, interaction: discord.Interaction, button: Button):
        user_key = (interaction.user.id, self.message_snowflake)
        now = datetime.now(timezone.utc)
        if user_key in self._last_called:
            delta = (now - self._last_called[user_key]).total_seconds()
            if delta < self.cooldown_seconds:
                reset_time = now + timedelta(seconds=(self.cooldown_seconds - delta))
                unix_timestamp = int(reset_time.timestamp())
                return await self.handler.send_message(interaction, content=f"Please wait <t:{unix_timestamp}:R> to fetch again.", ephemeral=True)
        self._last_called[user_key] = now
        state = cache.get(self.message_snowflake)
        if not state:
            return await self.handler.send_message(interaction, content="Timeout expired.")
        color = self.COLOR_MAP.get(state['health'].lower(), 0x5865F2)
        embed = discord.Embed(
            title="Information Statistics",
            color=color,
            timestamp=datetime.now(timezone.utc)
        )
        for key, value in state.get("kwargs", {}).items():
            embed.add_field(
                name=key,
                value=f"Valid: {value.get('valid')}\nPresent: {value.get('present')}",
                inline=False
            )
        embed.add_field(name="Date", value=state['date'].strftime("%Y-%m-%d %H:%M:%S UTC"), inline=True)
        embed.add_field(name="Executor", value=interaction.user.display_name, inline=True)
        embed.add_field(name="Speed", value=f"{state['speed']:.2f} sec.", inline=True)
        embed.add_field(name="Success", value=str(state["success"]), inline=True)
        await self.handler.send_message(interaction, embed=embed, ephemeral=True)

    @discord.ui.button(label="Report Issue", style=discord.ButtonStyle.secondary)
    async def report(self, interaction: discord.Interaction, button: Button):
        if interaction.user.id in self._reported_users:
            return await self.handler.send_message(interaction, content="You have already reported this message.")
        if interaction.user.id != interaction.user.id:
            return await self.handler.send_message(interaction, content="Only the command executor can report an issue.")
        self._reported_users.add(interaction.user.id)
        state = cache.get(self.message_snowflake)
        if not state:
            return await self.handler.send_message(interaction, content="Timeout expired.")
        try:
            message = await interaction.channel.fetch_message(self.message_snowflake)
            message_link = message.jump_url
        except (discord.NotFound, discord.Forbidden):
            logger.error("Failed to fetch message link.")
            return
        online_developers = []
        developers = await Developer.fetch_by_guild(guild_snowflake=interaction.guild.id)
        if developers:
            for dev in developers:
                member = interaction.guild.get_member(dev.member_snowflake)
                if member and member.status != discord.Status.offline:
                    online_developers.append(member.mention)
                try:
                    if member:
                        await member.send(f"Issue reported by {interaction.user.name}!\nMessage: {message_link}\n")
                except discord.Forbidden:
                    logger.warning(f"Unable to DM {member}.")
        if online_developers:
            msg = f"{', '.join(online_developers)} {'are' if len(online_developers) > 1 else 'is'} online and is now aware of the issue."
        else:
            msg = "All developers are offline, but your report was saved."
        await self.handler.send_message(interaction, content=msg)

class State:
    def __init__(self, ctx_or_interaction: Union[commands.Context, discord.Interaction], **kwargs):
        self.bot = DiscordBot.get_instance()
        self.counter = TimeToComplete()
        self.handler = MessageService(self.bot, self.bot.db_pool)
        self.kwargs = kwargs
        self.ctx_or_interaction = ctx_or_interaction
        if isinstance(ctx_or_interaction, commands.Context):
            self.start_time = ctx_or_interaction.message.created_at
        elif isinstance(ctx_or_interaction, discord.Interaction):
            self.start_time = discord.utils.utcnow()
        else:
            raise TypeError("Expected commands.Context or discord.Interaction")

    async def end(self, *, error: Union[str, discord.Embed, list[discord.Embed], discord.File] = None, success: Union[str, discord.Embed, list[discord.Embed], discord.File] = None, warning: Union[str, discord.Embed, list[discord.Embed], discord.File] = None, allowed_mentions: discord.AllowedMentions = discord.AllowedMentions.none()):
        message_obj = error or warning or success
        end_time = discord.utils.utcnow()
        elapsed = self.counter.time_elapsed_measurement(self.start_time, end_time)
        health = ""
        is_success = False
        if error:
            health = "\U0001F6AB"
        elif warning:
            health = "\U0001F950"
        elif success:
            health = "\U0001F3C6"
            is_success = True
        if isinstance(message_obj, list):
            if not message_obj:
                embed = discord.Embed(title="\u26A0\uFE0F", description="No results found.", color=discord.Color.red())
                return await self.handler.send_message(self.ctx_or_interaction, embed=embed, allowed_mentions=allowed_mentions)
            paginator = Paginator(self.handler.bot, self.ctx_or_interaction, message_obj)
            message = await paginator.start()
            cache[message.id] = {
                "date": self.start_time,
                "health": health,
                "kwargs": self.kwargs,
                "speed": elapsed,
                "success": success is not None
            }
            view = Information(message)
            merged_view = paginator.combine_views(view)
            await message.edit(view=merged_view)
            return message
        state_data = {"date": self.start_time, "health": health, "kwargs": self.kwargs, "speed": elapsed, "success": is_success}
        content = embed = file = None
        if isinstance(message_obj, str):
            content = message_obj
        elif isinstance(message_obj, discord.Embed):
            embed = message_obj
        elif isinstance(message_obj, discord.File):
            file = message_obj
        else:
            raise TypeError("Message must be str, discord.Embed, or discord.File")
        message = await self.handler.send_message(self.ctx_or_interaction, content=content, embed=embed, file=file, allowed_mentions=allowed_mentions)
        cache[message.id] = state_data
        view = Information(message)
        await message.edit(view=view)
        return message
