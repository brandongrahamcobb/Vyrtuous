
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
from contextlib import ExitStack
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

@pytest.fixture(scope="function")
def guild(bot):
    guild_obj = make_mock_guild(
        bot=bot,
        channel_defs={
            VOICE_CHANNEL_ONE_ID: (
                VOICE_CHANNEL_ONE_NAME,
                discord.ChannelType.voice
            ),
            VOICE_CHANNEL_TWO_ID: (
                VOICE_CHANNEL_TWO_NAME,
                discord.ChannelType.voice
            ),
            TEXT_CHANNEL_ID: (
                TEXT_CHANNEL_NAME,
                discord.ChannelType.text
            )
        },
        id=GUILD_ID,
        name=GUILD_NAME,
        members={},
        owner_id=PRIVILEGED_AUTHOR_ID,
        roles={
            ROLE_ID: SimpleNamespace(
                id=ROLE_ID,
                members=[],
                mention=f"<@&{ROLE_ID}>",
                name=ROLE_NAME
            )
        }
    )
    
    privileged_author_obj = make_mock_member(
        bot=True,
        guild=guild_obj,
        id=PRIVILEGED_AUTHOR_ID,
        name=PRIVILEGED_AUTHOR_NAME
    )
    
    not_privileged_author_obj = make_mock_member(
        bot=False,
        guild=guild_obj,
        id=NOT_PRIVILEGED_AUTHOR_ID,
        name=NOT_PRIVILEGED_AUTHOR_NAME
    )
    
    guild_obj.add_member(privileged_author_obj)
    guild_obj.add_member(not_privileged_author_obj)
    guild_obj.me = privileged_author_obj

    return guild_obj

@pytest_asyncio.fixture(scope="function")
async def bot():
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
    for cog in DISCORD_COGS:
        if cog != "vyrtuous.cogs.scheduled_tasks":
            await bot.load_extension(cog)
    type(bot).guilds = property(lambda self: [guild])
    bot._state = make_mock_state()
    yield bot
    await db_pool.close()

@pytest_asyncio.fixture(scope="function")
async def client():
    database: Optional[str] = os.getenv('POSTGRES_DB')
    host: Optional[str] = os.getenv('POSTGRES_HOST')
    password: Optional[str] = os.getenv('POSTGRES_PASSWORD')
    user: Optional[str] = os.getenv('POSTGRES_USER')
    dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"
    if not all([database, host, password, user]):
        pytest.skip("Database environment variables not set")
    db_pool = await asyncpg.create_pool(dsn=dsn)
    config = Config().get_config()
    client = DiscordClient(config=config, db_pool=db_pool)
    type(bot).guilds = property(lambda self: [guild_obj])
    yield client
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

async def make_capturing_send(channel, author):
    async def capturing_send(self, ctx, allowed_mentions=None, content=None, embed=None, file=None, **kwargs): 
        channel.messages.append({
            'content': content,
            'embed': embed,
            'file': file,
            'id': MESSAGE_ID
        })
        return make_mock_message(
            allowed_mentions=allowed_mentions,
            author=author,
            channel=channel,
            content=content,
            embeds=[embed],
            guild=ctx.guild,
            id=MESSAGE_ID
        )
    return capturing_send

async def edit(self, **kwargs):
    for k, v in kwargs.items():
        setattr(self, k, v)
    return self

async def prepared_command_handling(author, bot, channel, cog, content, guild, isinstance_patch, prefix):
    mock_message = make_mock_message(allowed_mentions=True, author=author, channel=channel, content=content, embeds=[], guild=guild, id=MESSAGE_ID)
    view = cmd_view.StringView(f"{prefix}{mock_message.content}")
    view.skip_string(prefix)
    mock_bot_user = guild.me
    with patch.object(type(bot), "user", new_callable=PropertyMock) as mock_user:
        mock_user.return_value = mock_bot_user
        view = cmd_view.StringView(mock_message.content)
        view.skip_string(prefix)
        ctx = Context(
            bot=bot,
            message=mock_message,
            prefix=prefix,
            view=view
        )
        command_name = content.lstrip(prefix).split()[0]
        ctx.command = bot.get_command(command_name)
        ctx.send = channel.send.__get__(channel, type(channel))
        ctx.invoked_with = command_name
        view.skip_string(command_name)
        view.skip_ws()
        fake_channels = {}
        if channel.type == discord.ChannelType.text:
            fake_channels = {
                guild.id: [
                    {"channel_id": channel.id, "channel_name": channel.name, "enabled": False}
                ]
            }
        cog_instance = bot.get_cog(cog)
        with ExitStack() as stack:
            if cog == "OwnerCommands":
                stack.enter_context(patch(
                    "vyrtuous.service.check_service.is_owner",
                    new=AsyncMock(return_value=True)
                ))
            capturing_send = await make_capturing_send(channel, author)
            cog_instance.handler.send_message = capturing_send.__get__(cog_instance.handler)
            def mock_isinstance(obj, cls):
                if cls == discord.VoiceChannel:
                    return hasattr(obj, 'type') and obj.type == discord.ChannelType.voice
                elif cls == discord.TextChannel:
                    return hasattr(obj, 'type') and obj.type == discord.ChannelType.text
                else:
                    return isinstance(obj, cls)
            stack.enter_context(patch.object(bot, "load_extension", new_callable=AsyncMock))
            stack.enter_context(patch.object(bot, "reload_extension", new_callable=AsyncMock))
            stack.enter_context(patch.object(bot, "unload_extension", new_callable=AsyncMock))
            if hasattr(cog_instance, "channel_service"):
                stack.enter_context(patch.object(cog_instance.channel_service, "resolve_channel", return_value=channel))
            if hasattr(cog_instance, "member_service"):
                def smart_resolve_member(ctx, member_id_or_mention):
                    if member_id_or_mention is None:
                        return None
                    raw = str(member_id_or_mention).strip()
                    if raw.startswith("<@") and raw.endswith(">"):
                        raw = raw.replace("<@", "").replace(">", "").lstrip("!")
                    if not raw.isdigit():
                        return None
                    return ctx.guild.get_member(int(raw))
                if hasattr(cog_instance, "member_service"):
                    stack.enter_context(patch.object(cog_instance.member_service, "resolve_member", side_effect=smart_resolve_member))
            stack.enter_context(patch(isinstance_patch, side_effect=mock_isinstance))
            stack.enter_context(patch(
                "vyrtuous.bot.discord_bot.DiscordBot.tree",
                new_callable=PropertyMock,
                return_value=MagicMock(
                    sync=AsyncMock(return_value=[]),
                    copy_global_to=MagicMock(),
                    clear_commands=MagicMock()
                )
            ))
            async def fetch_message(message_id):
                for msg in channel.messages:
                    if msg["id"] == message_id:
                        async def delete(self=None):
                            if self is None:
                                self = msg
                            channel.messages.remove(self)
                            return self
                        return SimpleNamespace(**msg, delete=delete)
            channel.fetch_message = fetch_message.__get__(channel)
            stack.enter_context(patch.object(channel, "fetch_message", new=fetch_message))
            stack.enter_context(patch("discord.utils.get", side_effect=lambda iterable, name=None: next((r for r in iterable if getattr(r, "name", None) == name), None )))
            await bot.invoke(ctx)