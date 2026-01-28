"""state.py A utility module for managing the failure, warning and success messages along with statistics about the bot's responses.

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

from typing import Union
from datetime import datetime, timezone
import asyncio
import uuid

from cachetools import TTLCache
from discord.ext import commands
import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.bug import Bug
from vyrtuous.db.roles.developer import Developer
from vyrtuous.service.message_service import MessageService, PaginatorService
from vyrtuous.utils.logger import logger
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.time_to_complete import TimeToComplete

cache = TTLCache(maxsize=500, ttl=8 * 60 * 60)


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
        source: Union[commands.Context, discord.Interaction, discord.Message],
    ):
        self.bot = DiscordBot.get_instance()
        self.config = self.bot.config
        self.message_service = MessageService(self.bot)
        self.counter = TimeToComplete()
        self.source = source
        self._reported_users = set()
        self.is_ephemeral = False
        self.start_time = self._get_start_time(source)
        self.message = None
        self.paginator = None
        if isinstance(source, commands.Context):
            self.submitter_id = source.author.id
        elif isinstance(source, discord.Interaction):
            self.submitter_id = source.user.id
        elif isinstance(source, discord.Message):
            self.submitter_id = source.author.id

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
            logger.warning(warning)
        elif error is not None:
            message_obj = error
            is_success = False
            logger.warning(error)
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
            self.paginator = PaginatorService(
                bot=self.bot,
                channel_source=self.source,
                pages=message_obj,
            )
            self.message = await self.paginator.start()
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
                content = f"{get_random_emoji()} {success}"
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
        if isinstance(self.source, discord.Interaction) and not paginated:
            if not self.source.response.is_done():
                await self.source.response.defer(ephemeral=True)
                self.is_ephemeral = True
            return await self.source.followup.send(**kwargs)
        elif isinstance(self.source, discord.Interaction) and paginated:
            if not self.source.response.is_done():
                await self.source.response.defer(ephemeral=True)
            return await self.source.followup.send(**kwargs)
        else:
            return await self.source.channel.send(**kwargs)

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
        self.bot.loop.create_task(self._wait_for_reactions())

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
                reaction, user = await self.bot.wait_for(
                    "reaction_add", timeout=30.0, check=look
                )
            except asyncio.TimeoutError:
                try:
                    await self.message.clear_reactions()
                except Exception as e:
                    logger.warning(str(e).capitalize())
                break
            await self._handle_reaction(reaction, user)
            try:
                await self.message.remove_reaction(reaction.emoji, user)
            except Exception as e:
                logger.warning(str(e).capitalize())

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
        if isinstance(self.source, discord.Interaction):
            await self.source.followup.send(
                f"{user.mention}, here is the info", embed=embed, ephemeral=True
            )
        else:
            try:
                await user.send(embed=embed)
            except discord.Forbidden as e:
                logger.info(str(e).capitalize())

    async def report_issue(self, user):
        reference = None
        if user.id in self._reported_users:
            try:
                await user.send("You already reported this message.")
            except discord.Forbidden as e:
                logger.info(str(e).capitalize())
            return
        self._reported_users.add(user.id)
        try:
            reference = str(uuid.uuid4())
            bug = Bug(
                channel_snowflake=self.message.channel.id,
                member_snowflakes=[],
                guild_snowflake=self.message.guild.id,
                id=reference,
                message_snowflake=self.message.id,
            )
            await bug.create()
        except discord.Forbidden as e:
            logger.info(str(e).capitalize())
        online_developer_mentions = []
        member = self.bot.get_user(self.config["discord_owner_id"])
        online_developer_mentions.append(member.mention)
        if self.source.guild:
            developers = await Developer.select(guild_snowflake=self.source.guild.id)
            message = f"Issue reported by {user.name}!\n**Message:** {self.message.jump_url}\n**Reference:** {reference}"
            for dev in developers:
                member = self.source.guild.get_member(dev.member_snowflake)
                if member and member.status != discord.Status.offline:
                    online_developer_mentions.append(member.mention)
                    try:
                        await member.send(message)
                    except discord.Forbidden as e:
                        logger.warning(
                            f"Unable to send a developer log ID: {id}. {str(e).capitalize()}"
                        )
        message = "Your report has been submitted"
        if online_developer_mentions:
            message = f"{message}. The developers {', '.join(online_developer_mentions)} are online and will respond to your report shortly."
        await user.send(message)

    async def send_pages(plural, pages, state):
        if pages:
            try:
                return await state.end(success=pages)
            except Exception as e:
                logger.info(str(e).capitalize())
                return await state.end(
                    warning="Embed size is too large. Limit the scope."
                )
        else:
            return await state.end(warning=f"No {plural} found.")
