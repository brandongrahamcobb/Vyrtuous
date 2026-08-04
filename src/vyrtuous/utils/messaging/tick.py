"""!/bin/python3
tick.py The purpose of this program is to deterministically report success, failure and warning metadata to messages.

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
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone

import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.cache.registry import MessageHistoryState
from vyrtuous.inc.helpers import resolve_author
from vyrtuous.upload import upload_service
from vyrtuous.utils.errors.error import CheckFailure
from vyrtuous.utils.messaging import emojis, message_service
from vyrtuous.utils.messaging.paginator import Paginator
from vyrtuous.utils.tracking import bug_service

COLOR_MAP = {"\u2705": 0x57F287, "\u26a0\ufe0f": 0xFEE65C, "\u274c": 0xED4245}

STATE_EMOJIS = {
    "\u2b05\ufe0f": -1,
    "\u27a1\ufe0f": 1,
    "\u2139\ufe0f": "info",
    "\U0001f4dd": "report",
    "\U0001f4c4": "upload",
}


@dataclass(frozen=True)
class TickRecord:
    date: datetime
    health_type: str
    speed: float
    success: bool


class Tick:
    def __init__(
        self,
        bot: DiscordBot,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        if ctx is None and interaction is None and message is None:
            raise CheckFailure("Discord source not defined to start message tick.")
        self.__bot = bot
        self.source = ctx or interaction or message
        self.author = resolve_author(self.source)
        self.elapsed = 0.0
        self.success = False
        self.start_time = self._resolve_start_time()

    @property
    def health_type(self) -> str:
        if self.elapsed <= 2.0:
            return "\u2705"
        elif self.elapsed <= 5.0:
            return "\u26a0\ufe0f"
        return "\u274c"

    def _resolve_start_time(self) -> datetime:
        if isinstance(self.source, commands.Context):
            return self.source.message.created_at
        elif isinstance(self.source, discord.Interaction):
            return discord.utils.utcnow()
        elif isinstance(self.source, discord.Message):
            return self.source.created_at
        raise TypeError("Expected Context, Interaction, or Message")

    def _compute_elapsed(self) -> float:
        return (discord.utils.utcnow() - self.start_time).total_seconds()

    def _resolve_content(
        self, success, warning, error
    ) -> tuple[str | discord.Embed | discord.File | list | None, bool]:
        if success is not None:
            return success, True
        elif warning is not None:
            self.__bot.logger.warning(warning)
            return warning, False
        elif error is not None:
            self.__bot.logger.warning(error)
            return error, False
        return None, True

    def _build_record(self, is_success: bool) -> TickRecord:
        return TickRecord(
            date=self.start_time,
            health_type=self.health_type,
            speed=self.elapsed,
            success=is_success,
        )

    async def _add_reactions(
        self, response: discord.Message, show_error_emoji: bool, paginated: bool
    ) -> None:
        if response.flags.ephemeral:
            return None
        if paginated:
            await response.add_reaction("\u2b05\ufe0f")
            await response.add_reaction("\u27a1\ufe0f")
        if not bool(self.__bot.config["release_mode"]):
            await response.add_reaction("\U0001f4c4")
        await response.add_reaction("\u2139\ufe0f")
        if show_error_emoji:
            await response.add_reaction("\U0001f4dd")
        self.__bot.loop.create_task(self._wait_for_reactions(response=response))

    async def _wait_for_reactions(self, response: discord.Message) -> None:
        def check(reaction, user):
            return (
                reaction.message.id == response.id
                and str(reaction.emoji) in STATE_EMOJIS
                and not user.bot
                and user.id == self.author.id
            )

        while True:
            try:
                reaction, user = await self.__bot.wait_for(
                    "reaction_add", timeout=30.0, check=check
                )
            except asyncio.TimeoutError:
                try:
                    await response.clear_reactions()
                except Exception as e:
                    self.__bot.logger.warning(str(e).capitalize())
                break
            await self._handle_reaction(response=response, reaction=reaction, user=user)
            try:
                await response.remove_reaction(reaction.emoji, user)
            except Exception as e:
                self.__bot.logger.warning(str(e).capitalize())

    async def _handle_reaction(
        self, response: discord.Message, reaction, user=None
    ) -> None:
        action = STATE_EMOJIS[str(reaction.emoji)]
        if action == "info":
            await self._show_info(user=user)
        elif action == "report":
            await self._report_issue(response=response, user=user)
        elif action == "upload":
            await upload_service.request_upload(source=self.source)

    async def _show_info(self, user) -> None:
        color = COLOR_MAP.get(self.health_type, 0xED4245)
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
        if isinstance(self.source, discord.Interaction):
            await self.source.followup.send(
                f"{user.mention}, here is the info", embed=embed, ephemeral=True
            )
        else:
            try:
                await user.send(embed=embed)
            except discord.Forbidden as e:
                self.__bot.logger.debug(str(e).capitalize())

    async def _report_issue(self, response: discord.Message, user) -> None:
        message_history = self.__bot.registry.get(MessageHistoryState)
        if user.id in message_history.reporters.get(response.id, []):
            try:
                await user.send("You already reported this message.")
            except discord.Forbidden as e:
                self.__bot.logger.debug(str(e).capitalize())
            return
        message_history.reporters[response.id].add(user.id)
        reference = str(uuid.uuid4())
        await bug_service.create_bug(message=response, reference=reference)
        # await developer_service.report_issue(
        #     author=self.author, message=response, reference=reference
        # )

    async def end(
        self,
        *,
        success=None,
        warning=None,
        error=None,
        allowed_mentions=discord.AllowedMentions.none(),
        view=None,
        ephemeral: bool = False,
    ) -> discord.Message:
        if self.source is None:
            raise CheckFailure("Source not provided.")
        message_obj, is_success = self._resolve_content(success, warning, error)
        self.elapsed = self._compute_elapsed()
        self.success = is_success
        show_error_emoji = bool(error or warning)

        registry = self.__bot.registry
        message_history = registry.get(MessageHistoryState)

        if isinstance(message_obj, list) and message_obj:
            paginator = Paginator()
            response = await paginator.start_by_message(
                ephemeral=ephemeral,
                pages=message_obj,
                source=self.source,
            )
            message_history.cache[response.id] = self._build_record(is_success)
            await self._add_reactions(
                response=response, show_error_emoji=show_error_emoji, paginated=True
            )
            return response

        content = embed = file = None
        if isinstance(message_obj, str):
            if success:
                content = f"{emojis.get_random_emoji()} {success}"
            elif warning:
                content = f"\u26a0\ufe0f {warning}"
            elif error:
                content = f"\u274c {error}"
        elif isinstance(message_obj, discord.Embed):
            embed = message_obj
        elif isinstance(message_obj, discord.File):
            file = message_obj
        elif message_obj is not None:
            raise TypeError("Message must be str, Embed, File, or list for pagination")
        response = await message_service.send_message(
            source=self.source,
            content=content,
            embed=embed,
            file=file,
            view=view,
            allowed_mentions=allowed_mentions,
            ephemeral=ephemeral,
        )
        message_history.cache[response.id] = self._build_record(is_success)
        await self._add_reactions(
            response=response, show_error_emoji=show_error_emoji, paginated=False
        )
        return response

    def update_source(
        self,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        self.source = ctx or interaction or message


# async def report_issue(author, message, reference) -> None:
#     bot: DiscordBot = DiscordBot.get_instance()
#     database_factory: DatabaseFactory = DatabaseFactory(MODEL)
#     online_developer_mentions = []
#     sysadmin = bot.get_user(bot.config["discord_owner_id"])
#     if sysadmin:
#         online_developer_mentions.append(sysadmin.mention)
#     if message.guild:
#         msg = f"Issue reported by {author.name}!\n**Message:** {message.jump_url}\n**Reference:** {reference}"
#         developers = await database_factory.select(singular=False)
#         for developer in developers:
#             member = message.guild.get_member(developer.member_snowflake)
#             if member and member.status != discord.Status.offline:
#                 online_developer_mentions.append(member.mention)
#                 try:
#                     await member.send(msg)
#                     if sysadmin:
#                         await sysadmin.send(msg)
#                 except discord.Forbidden as e:
#                     bot.logger.warning(
#                         f"Unable to send a developer log ID: {id}. {str(e).capitalize()}"
#                     )
#     msg = "Your report has been submitted"
#     if online_developer_mentions:
#         msg = f"{message}. The developers {', '.join(online_developer_mentions)} are online and will respond to your report shortly."
#     await author.send(msg)
