
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
from types import SimpleNamespace
from vyrtuous.inc.helpers import *
import asyncio

def make_mock_state():

    async def send_message(channel_id, content=None, embeds=None, **kwargs):
        return {
            'author': {
                'id': PRIVILEGED_AUTHOR_ID,
                'username': PRIVILEGED_AUTHOR_NAME
            },
            'channel_id': channel_id,
            'content': content,
            'embeds': embeds,
            'id': MESSAGE_ID,
        }

    mock_http = SimpleNamespace(
        allowed_mentions=None,
        send_message=send_message
    )

    return SimpleNamespace(
        allowed_mentions=None,
        http=mock_http,
        create_message=make_mock_message
    )

def make_mock_member(id, name, bot=True, voice_channel=False):

    async def edit(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        return self
    
    return type(
        'MockMember',
        (),
        {
            'bot': bot,
            'edit': edit,
            'id': id,
            'name': name,
            'mention': f'<@{id}>',
            'voice': SimpleNamespace(channel=voice_channel, mute=False)
        }
    )()

def make_mock_message(allowed_mentions, author, content, channel, embeds, guild, id):

    async def add_reaction(self, emoji):
        self.reactions.append(emoji)
    
    async def remove_reaction(self, emoji, user):
        if emoji in self.reactions:
            self.reactions.remove(emoji)
    
    async def clear_reactions(self):
        self.reactions.clear()
    
    async def edit(self, *, embed=None, content=None):
        if embed:
            self.edited_embeds.append(embed)
        if content is not None:
            self.content = content
        return self
    
    return type(
        'MockMessage',
        (),
        {
            'add_reaction': add_reaction,
            'allowed_mentions': allowed_mentions,
            'attachments': [],
            'author': author,
            'content': content,
            'channel': channel,
            'clear_reactions': clear_reactions,
            'edit': edit,
            'edited_embeds': [],
            'embeds': embeds,
            'guild': guild,
            'id': id,
            'reactions': [],
            'remove_reaction': remove_reaction,
            '_state': make_mock_state()
        }
    )()

def make_mock_guild(id, name, members, channel_defs, owner_id, roles):
    guild = type(
        'MockGuild',
        (),
        {
            'id': id,
            '_channels': {},
            'get_channel': lambda self, channel_id: self._channels.get(channel_id),
            'me': members.get(PRIVILEGED_AUTHOR_ID),
            '_members': members,
            'get_member': lambda self, member_id: self._members.get(member_id),
            'name': name,
            'owner_id': owner_id,
            "get_role": lambda self, role_id: roles.get(role_id),
        }
    )()
    channels = {}
    for cid, (name, channel_type) in channel_defs.items():
        channels[cid] = make_mock_channel(cid, name, guild, channel_type)
    guild._channels = channels
    for member in members.values():
        client_channel = list(channels.values())[0]
        member.voice.channel = client_channel
        client_channel.members.append(member)
    return guild

def make_mock_channel(id, name, guild, channel_type=None):

    async def async_send(self, content=None, embed=None, embeds=None, allowed_mentions=None, **kwargs):
        msg = make_mock_message(
            allowed_mentions=allowed_mentions,
            author=guild._members.get(list(guild._members.keys())[0]) if guild._members else None,
            content=content,
            embeds=embeds,
            guild=guild,
            id=id,
            _state=make_mock_state()
        )
        self.messages.append(msg)
        return msg

    return type(
        'MockChannel',
        (),
        {
            'guild': guild,
            'id': id,
            'members': [],
            'mention': f'<@{id}>',
            'messages': [],
            'name': name,
            'send': async_send,
            'type': channel_type,
            '_state': make_mock_state()
        }
    )()

async def mock_wait_for(event, timeout=None, check=None):
    raise asyncio.TimeoutError()