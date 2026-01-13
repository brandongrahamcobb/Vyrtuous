"""test_suite.py The purpose of this program is to provide the shared test variables for tests using Discord objects.

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
"""

import builtins
from contextlib import asynccontextmanager, contextmanager, ExitStack
from datetime import datetime, timezone
from discord.ext.commands import Context, view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import AsyncMock, MagicMock, PropertyMock, patch
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bot.discord_client import DiscordClient
from vyrtuous.config import Config
from vyrtuous.inc.helpers import *
from vyrtuous.database.roles.administrator import Administrator
from vyrtuous.database.roles.coordinator import Coordinator
from vyrtuous.database.roles.developer import Developer
from vyrtuous.database.roles.moderator import Moderator
from vyrtuous.utils.permission import PERMISSION_TYPES
from vyrtuous.service.state_service import State
from vyrtuous.tests.black_box.make_mock_objects import *
import asyncpg
import discord
import os
import pytest
import pytest_asyncio

RED = "\033[91m"
YELLOW = "\033[93m"
GREEN = "\033[92m"
RESET = "\033[0m"


@pytest.fixture(scope="function")
def guild():

    privileged_author_obj = create_member(
        bot=True, id=PRIVILEGED_AUTHOR_ID, name=PRIVILEGED_AUTHOR_NAME
    )

    not_privileged_author_obj = create_member(
        bot=False, id=NOT_PRIVILEGED_AUTHOR_ID, name=NOT_PRIVILEGED_AUTHOR_NAME
    )

    guild_obj = create_guild(
        bot=None,
        channels={},
        id=GUILD_ID,
        name=GUILD_NAME,
        members=[privileged_author_obj, not_privileged_author_obj],
        owner_snowflake=PRIVILEGED_AUTHOR_ID,
        roles={},
    )

    text_channel_obj = create_channel(
        channel_type="text",
        guild=guild_obj,
        id=TEXT_CHANNEL_ID,
        name=TEXT_CHANNEL_NAME,
        object_channel=discord.TextChannel,
    )

    voice_channel_one_obj = create_channel(
        channel_type="voice",
        guild=guild_obj,
        id=VOICE_CHANNEL_ONE_ID,
        name=VOICE_CHANNEL_ONE_NAME,
        object_channel=discord.VoiceChannel,
    )

    voice_channel_two_obj = create_channel(
        channel_type="voice",
        guild=guild_obj,
        id=VOICE_CHANNEL_TWO_ID,
        name=VOICE_CHANNEL_TWO_NAME,
        object_channel=discord.VoiceChannel,
    )

    channels = {
        TEXT_CHANNEL_ID: text_channel_obj,
        VOICE_CHANNEL_ONE_ID: voice_channel_one_obj,
        VOICE_CHANNEL_TWO_ID: voice_channel_two_obj,
    }

    guild_obj.channels = channels
    for member in guild_obj.members:
        voice_channel = list(guild_obj.channels.values())[1]
        member.voice.channel = voice_channel
        voice_channel.members.append(member)
    guild_obj.me = privileged_author_obj

    role_obj = create_role(
        guild=guild_obj, id=ROLE_ID, name=ROLE_NAME, members=guild_obj.members
    )
    guild_obj.roles = [role_obj]

    return guild_obj


@pytest_asyncio.fixture(scope="function")
async def bot(guild, privileged_author):
    database: Optional[str] = os.getenv("POSTGRES_DB")
    host: Optional[str] = os.getenv("POSTGRES_HOST")
    password: Optional[str] = os.getenv("POSTGRES_PASSWORD")
    user: Optional[str] = os.getenv("POSTGRES_USER")
    dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"
    if not all([database, host, password, user]):
        pytest.skip("Database environment variables not set")
    db_pool = await asyncpg.create_pool(dsn=dsn)
    config = Config().get_config()
    bot = DiscordBot(config=config, db_pool=db_pool)
    type(bot).guilds = property(lambda self: [guild])
    for cog in DISCORD_COGS:
        if cog != "vyrtuous.cogs.scheduled_tasks":
            await bot.load_extension(cog)
    yield bot
    await db_pool.close()


@pytest.fixture(scope="function")
def voice_channel_one(guild):
    return guild.channels[VOICE_CHANNEL_ONE_ID]


@pytest.fixture(scope="function")
def voice_channel_two(guild):
    return guild.channels[VOICE_CHANNEL_TWO_ID]


@pytest.fixture(scope="function")
def text_channel(guild):
    return guild.channels[TEXT_CHANNEL_ID]


@pytest.fixture(scope="function")
def role(guild):
    roles_list = guild.roles
    for role in roles_list:
        if role.id == ROLE_ID:
            return role


@pytest.fixture(scope="function")
def privileged_author(guild):
    return guild.get_member(PRIVILEGED_AUTHOR_ID)


@pytest.fixture(scope="function")
def not_privileged_author(guild):
    return guild.get_member(NOT_PRIVILEGED_AUTHOR_ID)


@pytest.fixture(scope="function")
def config(bot):
    config = bot.config
    yield config


@pytest.fixture(scope="function")
def prefix(config):
    prefix = config["discord_command_prefix"]
    yield prefix


async def prepare_context(bot, message, prefix):
    view = cmd_view.StringView(message.content)
    view.skip_string(prefix)
    command_name = view.get_word()
    ctx = Context(bot=bot, message=message, prefix=prefix, view=view)
    ctx.invoked_with = command_name
    ctx.command = bot.get_command(command_name)
    view.skip_ws()
    return ctx


def _normalize_payload(payload):
    if payload is None:
        return None, None, None
    if isinstance(payload, str):
        return payload, None, None
    if isinstance(payload, discord.Embed):
        return None, [payload], None
    if isinstance(payload, list):
        return None, payload, None
    if isinstance(payload, discord.File):
        return None, None, payload
    raise TypeError(f"Unsupported payload type: {type(payload)}")


@asynccontextmanager
async def capture(author, channel):
    captured = []
    send = State._send_message
    end = State.end

    async def _send(
        self,
        *,
        content=None,
        embed=None,
        embeds=None,
        file=None,
        paginated=False,
        allowed_mentions=None,
    ):
        if embed:
            embeds = [embed]
        return create_message(
            allowed_mentions=allowed_mentions,
            author=author,
            channel=channel,
            content=content,
            embeds=embeds,
            file=file,
            guild=channel.guild,
            id=MESSAGE_ID,
            paginated=paginated,
        )

    async def send(
        *,
        content=None,
        embed=None,
        embeds=None,
        file=None,
        paginated=False,
        allowed_mentions=None,
    ):
        if embed:
            embeds = [embed]
        return create_message(
            allowed_mentions=allowed_mentions,
            author=author,
            channel=channel,
            content=content,
            embeds=embeds,
            file=file,
            guild=channel.guild,
            id=MESSAGE_ID,
            paginated=paginated,
        )

    async def _end(self, *, success=None, warning=None, error=None, **kwargs):
        if success is not None:
            kind = "success"
            payload = success
        elif warning is not None:
            kind = "warning"
            payload = warning
        elif error is not None:
            kind = "error"
            payload = error
        else:
            kind = "unknown"
            content = None
        content, embeds, file = _normalize_payload(payload)
        embeds = [
            e.to_dict() if isinstance(e, discord.Embed) else e for e in (embeds or [])
        ]
        msg = await _send(self, content=content, embeds=embeds, file=file, **kwargs)
        captured.append({"type": kind, "message": msg})
        return msg

    # channel.send = _send
    State._send_message = _send
    channel.send = send
    State.end = _end
    try:
        yield captured
    finally:
        State.end = end
        State._send_message = send


@contextmanager
def prepare_discord_state(author, bot, channel, content, guild, highest_role):
    with ExitStack() as stack:
        stack.enter_context(
            patch(
                "discord.utils.get",
                side_effect=lambda iterable, name=None: next(
                    (r for r in iterable if r.name == name), None
                ),
            )
        )
        stack.enter_context(
            patch(
                "vyrtuous.service.check_service.has_equal_or_higher_role",
                new=AsyncMock(return_value=False),
            )
        )
        stack.enter_context(
            patch(
                "vyrtuous.service.check_service.moderator_predicator",
                new=AsyncMock(return_value=highest_role),
            )
        )
        # stack.enter_context(patch.object(bot, "get_channel", side_effect=lambda cid: channel if cid == channel.id else None))
        stack.enter_context(
            patch.object(
                bot,
                "get_guild",
                side_effect=lambda gid: guild if gid == guild.id else None,
            )
        )
        stack.enter_context(patch.object(bot, "load_extension", new_callable=AsyncMock))
        stack.enter_context(
            patch.object(bot, "reload_extension", new_callable=AsyncMock)
        )
        stack.enter_context(
            patch.object(bot, "unload_extension", new_callable=AsyncMock)
        )
        stack.enter_context(
            patch(
                "vyrtuous.service.paginator_service.Paginator.start",
                new_callable=AsyncMock,
            )
        )

        stack.enter_context(
            patch(
                "vyrtuous.bot.discord_bot.DiscordBot.tree",
                new_callable=PropertyMock,
                return_value=MagicMock(
                    sync=AsyncMock(return_value=[]),
                    copy_global_to=MagicMock(),
                    clear_commands=MagicMock(),
                ),
            )
        )
        yield


async def command(bot, ctx):
    bot.loop = asyncio.get_running_loop()
    await bot.invoke(ctx)


async def on_message(bot, message):
    bot.loop = asyncio.get_running_loop()
    bot.dispatch("message", message)
    await asyncio.sleep(1)


async def prepared_command_handling(
    author, bot, channel, content, guild, highest_role, prefix
):
    message = create_message(
        allowed_mentions=True,
        author=author,
        bot=bot,
        channel=channel,
        content=content,
        guild=guild,
        id=MESSAGE_ID,
    )
    mock_bot_user = guild.me
    with patch.object(
        bot, "_connection", create=True
    ) as mock_conn, ExitStack() as stack, prepare_discord_state(
        author, bot, channel, content, guild, highest_role
    ):
        mock_conn.user = mock_bot_user
        mock_conn.return_value = mock_bot_user
        ctx = await prepare_context(bot, message, prefix)
        async with capture(author, channel) as captured:
            if ctx.command:
                await command(bot, ctx)
            else:
                message.content = f"{prefix}{message.content}"
                await on_message(bot, message)
    return captured


def extract_embed_text(embed: discord.Embed) -> str:
    parts = []
    if embed.title:
        parts.append(embed.title)
    if embed.description:
        parts.append(embed.description)
    for field in embed.fields:
        parts.append(f"{field.name}: {field.value}")
    return "\n".join(parts)


@pytest_asyncio.fixture(scope="function")
async def permission(request, voice_channel_one, guild, privileged_author):
    perm_type = request.param
    if perm_type is None:
        yield None
        return
    PERMISSION_CLASSES = {
        "Administrator": Administrator,
        "Coordinator": Coordinator,
        "Developer": Developer,
        "Moderator": Moderator,
    }
    perm_class = PERMISSION_CLASSES[perm_type]
    CHANNEL_CLASSES = {"Coordinator", "Moderator"}
    ROLE_CLASSES = {"Administrator"}
    if perm_type in CHANNEL_CLASSES:
        perm_instance = perm_class(
            channel_snowflake=voice_channel_one.id,
            guild_snowflake=guild.id,
            member_snowflake=privileged_author.id,
        )
    elif perm_type in ROLE_CLASSES:
        perm_instance = perm_class(
            guild_snowflake=guild.id,
            member_snowflake=privileged_author.id,
            role_snowflakes=[ROLE_ID],
        )
    else:
        perm_instance = perm_class(
            guild_snowflake=guild.id, member_snowflake=privileged_author.id
        )
    await perm_instance.grant()
    try:
        yield perm_instance
    finally:
        await perm_instance.revoke()
