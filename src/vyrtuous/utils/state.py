
''' state.py A utility module for managing the failure, warning and success messages along with statistics about the bot's responses.

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
from typing import Union
from cachetools import TTLCache
from datetime import datetime, timezone
from discord.ext import commands
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.message_service import MessageService
from vyrtuous.utils.developer import Developer
from vyrtuous.utils.developer_log import DeveloperLog
from vyrtuous.utils.paginator import Paginator
from vyrtuous.utils.time_to_complete import TimeToComplete
import asyncio
import discord

cache = TTLCache(maxsize=500, ttl=8*60*60)

class State:

    COLOR_MAP = {
        "\u2705": 0x57F287,
        "\u26A0\ufe0f": 0xFEE65C,
        "\u274C": 0xED4245
    }

    STATE_EMOJIS = {
        "\u2b05\ufe0f": -1,
        "\u27a1\ufe0f": 1,
        "\u2139\ufe0f": "info",
        "\U0001F4DD": "report"
    }

    def __init__(self, ctx_or_interaction: Union[commands.Context, discord.Interaction, discord.Message]):
        self.bot = DiscordBot.get_instance()
        self.config = self.bot.config
        self.message_service = MessageService(self.bot, self.bot.db_pool)
        self.counter = TimeToComplete()
        self.ctx_or_interaction = ctx_or_interaction
        self._reported_users = set()
        self.is_ephemeral = False
        self.start_time = self._get_start_time(ctx_or_interaction)
        self.message = None
        self.paginator = None
        if isinstance(ctx_or_interaction, commands.Context):
            self.submitter_id = ctx_or_interaction.author.id
        elif isinstance(ctx_or_interaction, discord.Interaction):
            self.submitter_id = ctx_or_interaction.user.id
        elif isinstance(ctx_or_interaction, discord.Message):
            self.submitter_id = ctx_or_interaction.author.id

    def _get_start_time(self, ctx_or_interaction):
        if isinstance(ctx_or_interaction, commands.Context):
            return ctx_or_interaction.message.created_at
        elif isinstance(ctx_or_interaction, discord.Interaction):
            return discord.utils.utcnow()
        elif isinstance(ctx_or_interaction, discord.Message):
            return ctx_or_interaction.created_at
        else:
            raise TypeError("Expected Context, Interaction, or Message")

    async def end(self, *, success=None, warning=None, error=None, message_obj=None, allowed_mentions=discord.AllowedMentions.none()):
        if success is not None:
            message_obj = success
            is_success = True
        elif warning is not None:
            message_obj = warning
            is_success = False
        elif error is not None:
            message_obj = error
            is_success = False
        else:
            message_obj = None
            is_success = True
        elapsed = self.counter.time_elapsed_measurement(self.start_time, discord.utils.utcnow())
        if elapsed <= 2.0:
            health_type = "\u2705"
        elif elapsed <= 5.0:
            health_type = "\u26A0\ufe0f"
        else:
            health_type = "\u274C"
        self.health_type = health_type
        self.success = is_success
        self.elapsed = elapsed
        if isinstance(message_obj, list) and message_obj:
            self.paginator = Paginator(bot=self.bot, ctx_or_interaction=self.ctx_or_interaction, pages=message_obj)
            self.message = await self.paginator.start()
            cache[self.message.id] = {
                "date": self.start_time,
                "health_type": health_type,
                "speed": elapsed,
                "success": is_success
            }
            await self._add_reactions(paginated=True)
            return self.message
        content = embed = file = None
        if isinstance(message_obj, str):
            content = message_obj
        elif isinstance(message_obj, discord.Embed):
            embed = message_obj
        elif isinstance(message_obj, discord.File):
            file = message_obj
        else:
            raise TypeError("Message must be str, embed, file, or list for pagination")
        self.message = await self._send_message(content=content, embed=embed, file=file, paginated=False, allowed_mentions=allowed_mentions)
        cache[self.message.id] = {
            "date": self.start_time,
            "health_type": health_type,
            "speed": elapsed,
            "success": is_success
        }
        await self._add_reactions(paginated=False)
        return self.message

    async def _send_message(self, content=None, embed=None, file=None, paginated=False, allowed_mentions=discord.AllowedMentions.none()):
        kwargs = {
            "content": content,
            "embed": embed,
            "allowed_mentions": allowed_mentions,
        }
        if file is not None:
            kwargs["file"] = file
        if isinstance(self.ctx_or_interaction, discord.Interaction) and not paginated:
            if not self.ctx_or_interaction.response.is_done():
                await self.ctx_or_interaction.response.defer(ephemeral=True)
                self.is_ephemeral = True
            return await self.ctx_or_interaction.followup.send(**kwargs)
        elif isinstance(self.ctx_or_interaction, discord.Interaction) and paginated:
            if not self.ctx_or_interaction.response.is_done():
                await self.ctx_or_interaction.response.defer()
            return await self.ctx_or_interaction.followup.send(**kwargs)
        else:
            return await self.ctx_or_interaction.channel.send(**kwargs)


    async def _add_reactions(self, paginated: bool):
        if not self.message or self.is_ephemeral:
            return
        if self.message and self.message.webhook_id is not None:
            return
        if paginated:
            await self.message.add_reaction("\u2b05\ufe0f")
            await self.message.add_reaction("\u27a1\ufe0f")
        await self.message.add_reaction("\u2139\ufe0f")
        await self.message.add_reaction("\U0001F4DD")
        self.bot.loop.create_task(self._wait_for_reactions())

    async def _wait_for_reactions(self):
        def check(reaction, user):
            return (
                reaction.message.id == self.message.id
                and str(reaction.emoji) in self.STATE_EMOJIS
                and not user.bot
                and user.id == self.submitter_id
            )
        while True:
            try:
                reaction, user = await self.bot.wait_for("reaction_add", timeout=30.0, check=check)
            except asyncio.TimeoutError:
                try:
                    await self.message.clear_reactions()
                except:
                    pass
                break
            await self._handle_reaction(reaction, user)
            try:
                await self.message.remove_reaction(reaction.emoji, user)
            except:
                pass

    async def _handle_reaction(self, reaction, user=None):
        action = self.STATE_EMOJIS[str(reaction.emoji)]
        if action == "info":
            await self.show_info(user)
        elif action == "report":
            await self.report_issue(user)

    async def show_info(self, user):
        color = self.COLOR_MAP.get(self.health_type)
        embed = discord.Embed(title="Information Statistics", color=color, timestamp=datetime.now(timezone.utc))
        embed.add_field(name="Date", value=self.start_time.strftime("%Y-%m-%d %H:%M:%S UTC"), inline=True)
        embed.add_field(name="Health", value=self.health_type, inline=True)
        embed.add_field(name="Speed", value=f"{self.elapsed:.2f} sec.", inline=True)
        embed.add_field(name="Success", value=str(self.success), inline=True)
        if isinstance(self.ctx_or_interaction, discord.Interaction):
            await self.ctx_or_interaction.followup.send(f"{user.mention}, here is the info", embed=embed, ephemeral=True)
        else:
            try:
                await user.send(embed=embed)
            except discord.Forbidden:
                pass

    async def report_issue(self, user):
        reference = None
        if user.id in self._reported_users:
            try:
                await user.send("You already reported this message.")
            except discord.Forbidden:
                pass
            return
        self._reported_users.add(user.id)
        try:
            developer_log = DeveloperLog(channel_snowflake=self.message.channel.id, developer_snowflakes=[], guild_snowflake=self.message.guild.id, message_snowflake=self.message.id)
            reference = await developer_log.create()
            id = reference.id
            await user.send("Your report has been submitted.")
        except discord.Forbidden:
            pass
        if self.ctx_or_interaction.guild:
            developers = await Developer.fetch_by_guild(self.ctx_or_interaction.guild.id)
            for dev in developers:
                member = self.ctx_or_interaction.guild.get_member(dev.member_snowflake)
                if member and member.status != discord.Status.offline:
                    try:
                        await member.send(f"Issue reported by {user.name}!\nMessage: {self.message.jump_url}\n**Reference:** {id}")
                    except discord.Forbidden:
                        pass
                member = self.bot.get_user(self.config['discord_owner_id'])
                try:
                    await member.send(f"Issue reported by {user.name}!\n**Message:** {self.message.jump_url}\n**Reference:** {id}")
                except discord.Forbidden:
                    pass
