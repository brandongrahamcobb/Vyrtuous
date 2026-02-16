"""!/bin/python3
state_service.py The purpose of this program is to extend Service to service the success, failure and warning messages.

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

import asyncio
import uuid
from datetime import datetime, timezone

import discord
from cachetools import TTLCache
from discord.ext import commands

from vyrtuous.utils.message_service import PaginatorService

cache = TTLCache(maxsize=500, ttl=8 * 60 * 60)


class TimeToComplete:
    def is_around_one_second(self, elapsed: float = 1.0):
        return 0.0 <= elapsed <= 2.0

    def time_elapsed_measurement(self, start: datetime, end: datetime) -> float:
        if start is None or end is None:
            return 0.0
        return (end - start).total_seconds()


class StateService:
    COLOR_MAP = {"\u2705": 0x57F287, "\u26a0\ufe0f": 0xFEE65C, "\u274c": 0xED4245}

    STATE_EMOJIS = {
        "\u2b05\ufe0f": -1,
        "\u27a1\ufe0f": 1,
        "\u2139\ufe0f": "info",
        "\U0001f4dd": "report",
    }

    def __init__(
        self,
        *,
        author_service=None,
        bot=None,
        bug_service=None,
        developer_service=None,
        emoji=None,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        if (ctx is None) == (interaction is None) == (message is None):
            raise commands.CheckFailure("Discord source not defined in StateService.")
        self.__author_service = author_service
        self.__bot = bot
        self.__config = self.__bot.config
        self.__bug_service = bug_service
        self.__developer_service = developer_service
        self.__emoji = emoji
        self.counter = TimeToComplete()
        self._source = ctx or interaction or message
        self._reported_users = set()
        self.is_ephemeral = False
        self.start_time = self._get_start_time(self._source)
        self.message = discord.Message | None
        self.__paginator = None
        self.submitter_id = self.__author_service.resolve_author(self._source).id

    def _get_start_time(self, source):
        if isinstance(source, commands.Context):
            return source.message.created_at
        elif isinstance(source, discord.Interaction):
            return discord.utils.utcnow()
        elif isinstance(source, discord.Message):
            return source.created_at
        else:
            raise TypeError("Expected Context, Interaction, or Message")

    async def end(
        self,
        *,
        success=None,
        warning=None,
        error=None,
        message_obj=None,
        allowed_mentions=discord.AllowedMentions.none(),
    ):
        if success is not None:
            message_obj = success
            is_success = True
        elif warning is not None:
            message_obj = warning
            is_success = False
            self.__bot.logger.warning(warning)
        elif error is not None:
            message_obj = error
            is_success = False
            self.__bot.logger.warning(error)
        else:
            message_obj = None
            is_success = True
        elapsed = self.counter.time_elapsed_measurement(
            self.start_time, discord.utils.utcnow()
        )
        if elapsed <= 2.0:
            health_type = "\u2705"
        elif elapsed <= 5.0:
            health_type = "\u26a0\ufe0f"
        else:
            health_type = "\u274c"
        self.health_type = health_type
        self.success = is_success
        self.elapsed = elapsed
        show_error_emoji = False
        if error or warning:
            show_error_emoji = True
        if isinstance(message_obj, list) and message_obj:
            self.__paginator = PaginatorService(
                bot=self.__bot,
                channel_source=self._source,
                pages=message_obj,
            )
            self.message = await self.__paginator.start()
            cache[self.message.id] = {
                "date": self.start_time,
                "health_type": health_type,
                "speed": elapsed,
                "success": is_success,
            }
            await self._add_reactions(show_error_emoji=show_error_emoji, paginated=True)
            return self.message
        content = embed = file = None
        if isinstance(message_obj, str):
            if success:
                content = f"{self.__emoji.get_random_emoji()} {success}"
            if warning:
                content = f"\u26a0\ufe0f {warning}"
            if error:
                content = f"\u274c {error}"
        elif isinstance(message_obj, discord.Embed):
            embed = message_obj
        elif isinstance(message_obj, discord.File):
            file = message_obj
        else:
            raise TypeError("Message must be str, embed, file, or list for pagination")
        self.message = await self._send_message(
            content=content,
            embed=embed,
            file=file,
            paginated=False,
            allowed_mentions=allowed_mentions,
        )
        cache[self.message.id] = {
            "date": self.start_time,
            "health_type": health_type,
            "speed": elapsed,
            "success": is_success,
        }
        await self._add_reactions(show_error_emoji=show_error_emoji, paginated=False)
        return self.message

    async def _send_message(
        self,
        content=None,
        embed=None,
        file=None,
        paginated=False,
        allowed_mentions=discord.AllowedMentions.none(),
    ):
        kwargs = {
            "content": content,
            "embed": embed,
            "allowed_mentions": allowed_mentions,
        }
        if file is not None:
            kwargs["file"] = file
        if isinstance(self._source, discord.Interaction) and not paginated:
            if not self._source.response.is_done():
                await self._source.response.defer(ephemeral=True)
                self.is_ephemeral = True
            return await self._source.followup.send(**kwargs)
        elif isinstance(self._source, discord.Interaction) and paginated:
            if not self._source.response.is_done():
                await self._source.response.defer(ephemeral=True)
            return await self._source.followup.send(**kwargs)
        else:
            return await self._source.channel.send(**kwargs)

    async def _add_reactions(self, show_error_emoji: bool, paginated: bool):
        if not self.message or self.is_ephemeral:
            return
        if self.message and self.message.webhook_id is not None:
            return
        if paginated:
            await self.message.add_reaction("\u2b05\ufe0f")
            await self.message.add_reaction("\u27a1\ufe0f")
        await self.message.add_reaction("\u2139\ufe0f")
        if show_error_emoji:
            await self.message.add_reaction("\U0001f4dd")
        self.__bot.loop.create_task(self._wait_for_reactions())

    async def _wait_for_reactions(self):
        def look(reaction, user):
            return (
                reaction.message.id == self.message.id
                and str(reaction.emoji) in self.STATE_EMOJIS
                and not user.bot
                and user.id == self.submitter_id
            )

        while True:
            try:
                reaction, user = await self.__bot.wait_for(
                    "reaction_add", timeout=30.0, check=look
                )
            except asyncio.TimeoutError:
                try:
                    await self.message.clear_reactions()
                except Exception as e:
                    self.__bot.logger.warning(str(e).capitalize())
                break
            await self._handle_reaction(reaction, user)
            try:
                await self.message.remove_reaction(reaction.emoji, user)
            except Exception as e:
                self.__bot.logger.warning(str(e).capitalize())

    async def _handle_reaction(self, reaction, user=None):
        action = self.STATE_EMOJIS[str(reaction.emoji)]
        if action == "info":
            await self.show_info(user)
        elif action == "report":
            await self.report_issue(user)

    async def show_info(self, user):
        color = self.COLOR_MAP.get(self.health_type, "0xED4245")
        embed = discord.Embed(
            title="Information Statistics",
            color=color,
            timestamp=datetime.now(timezone.utc),
        )
        embed.add_field(
            name="Date",
            value=self.start_time.strftime("%Y-%m-%d %H:%M:%S UTC"),
            inline=True,
        )
        embed.add_field(name="Health", value=self.health_type, inline=True)
        embed.add_field(name="Speed", value=f"{self.elapsed:.2f} sec.", inline=True)
        embed.add_field(name="Success", value=str(self.success), inline=True)
        if isinstance(self._source, discord.Interaction):
            await self._source.followup.send(
                f"{user.mention}, here is the info", embed=embed, ephemeral=True
            )
        else:
            try:
                await user.send(embed=embed)
            except discord.Forbidden as e:
                self.__bot.logger.info(str(e).capitalize())

    async def report_issue(self, user):
        reference = None
        if user.id in self._reported_users:
            try:
                await user.send("You already reported this message.")
            except discord.Forbidden as e:
                self.__bot.logger.info(str(e).capitalize())
            return
        self._reported_users.add(user.id)
        reference = str(uuid.uuid4())
        await self.__bug_service.create_bug(message=self.message, reference=reference)
        await self.__developer_service.report_issue(
            message=self.message, reference=reference, source=self._source, user=user
        )

    @classmethod
    async def send_pages(cls, title, pages, state):
        if pages:
            return await state.end(success=pages)
        else:
            return await state.end(success=f"No {title.lower()} found.")
