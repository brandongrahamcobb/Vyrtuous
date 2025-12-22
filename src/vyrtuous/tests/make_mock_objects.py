
''' make_mock_objects.py The purpose of this program is to provide the shared test variables for tests using Discord objects.
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
import discord

def make_mock_state():

    async def send_message(author=None, channel_id=None, content=None, embeds=None, params=None, **kwargs):
        return {
            'author': {
                'id': author.id,
                'username': author.name
            },
            'channel_id': channel_id,
            'content': content,
            'embeds': embeds,
            'id': MESSAGE_ID,
            'params': params
        }

    async def create_message(channel=None, embeds=None, **kwargs):
        msg = make_mock_message(channel=channel, embeds=embeds)
        if hasattr(channel, "messages"):
            channel.messages.append(msg)
        return msg

    mock_http = SimpleNamespace(
        allowed_mentions=None,
        create_message=create_message,
        send_message=send_message
    )

    return SimpleNamespace(
        allowed_mentions=None,
        http=mock_http,
        create_message=make_mock_message
    )

def make_mock_member(bot=True, guild=None, id=None, name=None, voice_channel=False):

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
            'guild': guild,
            'id': id,
            'name': name,
            'mention': f'<@{id}>',
            'voice': SimpleNamespace(channel=voice_channel, mute=False)
        }
    )()

def make_mock_message(allowed_mentions=None, author=None, content=None, channel=None, embeds=None, guild=None, id=None, **kwargs):

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
            'embeds': embeds or [],
            'guild': guild,
            'id': id,
            'reactions': [],
            'remove_reaction': remove_reaction,
            '_state': make_mock_state()
        }
    )()

def make_mock_guild(bot, channel_defs=None, id=None, name=None, members=None, owner_id=None, roles=None):
    guild = type(
        'MockGuild',
        (),
        {
            'id': id,
            'channels': {},
            'get_channel': lambda self, channel_id: self.channels.get(channel_id),
            'me': bot,
            'members': members,
            'get_member': lambda self, member_id: self.members.get(member_id),
            'name': name,
            'owner_id': owner_id,
            "get_role": lambda self, role_id: roles.get(role_id),
        }
    )()
    channels = {}
    for cid, (name, channel_type) in channel_defs.items():
        channels[cid] = make_mock_channel(channel_type=channel_type, guild=guild, id=cid, name=name)
    guild.channels = channels
    for member in (members or {}).values():
        client_channel = list(channels.values())[0]
        member.voice.channel = client_channel
        client_channel.members.append(member)
        for role in roles.values():
            role.members.append(member)
    guild.add_member = lambda member: guild.members.update({member.id: member})
    return guild

def make_mock_channel(channel_type=None, guild=None, id=None, name=None):

    async def async_send(self, allowed_mentions=None, author=None, content=None, embed=None, embeds=None, file=None, **kwargs):
        self.messages.append({
            'content': content,
            'embed': embed,
            'file': file,
            'id': MESSAGE_ID
        })
        msg = make_mock_message(
            allowed_mentions=allowed_mentions,
            author=author,
            channel=self,
            content=content,
            embeds=embeds or ([embed] if embed else []),
            guild=guild,
            id=id
        )
        return msg

    return type(
        'MockChannel',
        (),
        {
            '_state': make_mock_state(),
            'guild': guild,
            'id': id,
            'members': [],
            'mention': f'<#{id}>',
            'messages': [],
            'name': name,
            'send': async_send,
            'send_messages': True,
            'type': channel_type
        }
    )()

async def mock_wait_for(event, timeout=None, check=None):
    raise asyncio.TimeoutError()