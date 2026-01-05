
''' test_suite.py The purpose of this program is to provide the shared test variables for tests using Discord objects.

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
import builtins
from contextlib import contextmanager, ExitStack
from datetime import datetime, timezone
from discord.ext.commands import Context, view as cmd_view
from types import SimpleNamespace
from typing import Optional
from unittest.mock import AsyncMock, MagicMock, PropertyMock, patch
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.bot.discord_client import DiscordClient
from vyrtuous.config import Config
from vyrtuous.inc.helpers import *
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
    text_channel_obj = create_channel(
        channel_type='text',
        id=TEXT_CHANNEL_ID,
        name=TEXT_CHANNEL_NAME,
        object_channel=discord.TextChannel
    )

    voice_channel_obj = create_channel(
        channel_type='voice'
        id=VOICE_CHANNEL_ONE_ID,
        name=VOICE_CHANNEL_ONE_NAME,
        object_channel=discord.VoiceChannel
    )

    channels = {TEXT_CHANNEL_ID: text_channel, VOICE_CHANNEL_ONE_ID: voice_channel}

    guild_obj = create_guild(
        bot=None,
        channels=channels,
        id=GUILD_ID,
        name=GUILD_NAME,
        members={},
        owner_id=PRIVILEGED_AUTHOR_ID,
        roles={}
    )

    for channel in guild_obj.channels.values():
        channel.guild = guild_obj
    
    privileged_author_obj = create_member(
        bot=True,
        guild=guild_obj,
        id=PRIVILEGED_AUTHOR_ID,
        name=PRIVILEGED_AUTHOR_NAME
    )
    
    not_privileged_author_obj = create_member(
        bot=False,
        guild=guild_obj,
        id=NOT_PRIVILEGED_AUTHOR_ID,
        name=NOT_PRIVILEGED_AUTHOR_NAME
    )
    
    guild_obj.members[privileged_author_obj.id] = privileged_author_obj
    guild_obj.members[not_privileged_author_obj.id] = not_privileged_author_obj
    for member in guild_obj.members.values():
        voice_channel = list(guild_obj.channels.values())[1]
        member.voice.channel = voice_channel
        voice_channel.members.append(member)
    guild_obj.me = privileged_author_obj
    
    role_obj = create_role(
        guild=guild_obj,
        id=ROLE_ID,
        name=ROLE_NAME,
        members=list(guild_obj.members.values()))
    )
    guild_obj.roles[ROLE_ID] = role_obj

    return guild_obj

@pytest_asyncio.fixture(scope="function")
async def bot(guild):
    database: Optional[str] = os.getenv('POSTGRES_DB')
    host: Optional[str] = os.getenv('POSTGRES_HOST')
    password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
    user: Optional[str] = os.getenv('POSTGRES_USER')
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
    roles_list = list(guild.roles.values())
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
    prefix = config['discord_command_prefix']
    yield prefix

async def prepare_context(bot, message, prefix):
    view = cmd_view.StringView(message.content)
    view.skip_string(prefix)
    command_name = view.get_word()
    ctx = Context(
        bot=bot,
        message=message,
        prefix=prefix,
        view=view
    )
    ctx.invoked_with = command_name
    ctx.command = bot.get_command(command_name)
    view.skip_ws()
    return ctx

async def capture(channel):
    captured = []
    send = State._send_message
    async def _send(self, content=None, embed=None, file=None, paginated=False, allowed_mentions=None):
        msg = await send(self, content=content, embed=embed, file=file, paginated=paginated, allowed_mentions=allowed_mentions)
        captured.append(msg)
        return msg
    State._send_message = _send
    try:
        yield captured
    finally:
        State._send_message = send

@contextmanager
def prepare_discord_state(bot, channel, guild, highest_role):
    with ExitStack() as stack:
        stack.enter_context(patch("discord.utils.get", side_effect=lambda iterable, name=None: next((r for r in iterable if r.name == name), None)))
        stack.enter_context(patch("vyrtuous.service.check_service.has_equal_or_higher_role", new=AsyncMock(return_value=False)))
        stack.enter_context(patch("vyrtuous.service.check_service.is_system_owner_developer_guild_owner_administrator_coordinator", new=AsyncMock(return_value=highest_role)))
        stack.enter_context(patch("vyrtuous.service.check_service.is_system_owner_developer_guild_owner_administrator_coordinator_moderator", new=AsyncMock(return_value=highest_role)))
        stack.enter_context(patch.object(bot, "get_channel", side_effect=lambda cid: guild.get_channel(cid)))
        stack.enter_context(patch.object(bot, "get_guild", side_effect=lambda gid: guild if gid == guild.id else None))
        stack.enter_context(patch.object(bot, "get_role", side_effect=lambda uid: guild.get_role(uid)))
        stack.enter_context(patch.object(bot, "get_user", side_effect=lambda uid: guild.get_member(uid)))
        stack.enter_context(patch.object(bot, "load_extension", new_callable=AsyncMock))
        stack.enter_context(patch.object(bot, "reload_extension", new_callable=AsyncMock))
        stack.enter_context(patch.object(bot, "unload_extension", new_callable=AsyncMock))
        yield

async def command(bot, ctx):
    await bot.invoke(ctx)

async def on_message(bot, message):
    bot.loop = asyncio.get_running_loop()
    bot.dispatch("message", message)
    await asyncio.sleep(0)

async def prepared_command_handling(author, bot, channel, content, guild, highest_role, prefix):
    captured = []
    message = create_message(
        allowed_mentions=True,
        author=author,
        channel=channel,
        content=content,
        guild=guild,
        id=MESSAGE_ID
    )
    mock_bot_user = guild.me
    with patch.object(bot, "_connection", create=True) as mock_conn, ExitStack() as stack:
        mock_conn.user = mock_bot_user
        mock_conn.return_value = mock_bot_user
        ctx = await prepare_context(bot, message, prefix)
        with capture(channel) as captured_messages, prepare_discord_state(bot, channel, guild, highest_role):
            await command(bot, ctx)
            await on_message(bot, message)
    return captured_messages
        
def extract_embed_text(embed: discord.Embed) -> str:
    parts = []
    if embed.title:
        parts.append(embed.title)
    if embed.description:
        parts.append(embed.description)
    for field in embed.fields:
        parts.append(f"{field.name}: {field.value}")
    return "\n".join(parts)

# def _get_start_time(self, ctx_or_interaction):
#     if hasattr(ctx_or_interaction, "created_at"):
#         return ctx_or_interaction.created_at
#     if hasattr(ctx_or_interaction, "message") and hasattr(ctx_or_interaction.message, "created_at"):
#         return ctx_or_interaction.message.created_at
#     return datetime.now(timezone.utc)

# @pytest_asyncio.fixture(scope="function")
# async def client():
#     database: Optional[str] = os.getenv('POSTGRES_DB')
#     host: Optional[str] = os.getenv('POSTGRES_HOST')
#     password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
#     user: Optional[str] = os.getenv('POSTGRES_USER')
#     dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"
#     if not all([database, host, password, user]):
#         pytest.skip("Database environment variables not set")
#     db_pool = await asyncpg.create_pool(dsn=dsn)
#     config = Config().get_config()
#     client = DiscordClient(config=config, db_pool=db_pool)
#     type(bot).guilds = property(lambda self: [guild])
#     yield client
#     await db_pool.close()