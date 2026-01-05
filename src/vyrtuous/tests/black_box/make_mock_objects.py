
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
from datetime import datetime, timezone
from types import SimpleNamespace
from vyrtuous.inc.helpers import *
import asyncio
import discord

def create_http():
    http = SimpleNamespace(allowed_mentions=None)
    return http

def create_member(bot=False, channel=None, guild=None, id=None, name=None):

    async def edit(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        return self

    async def update_voice_state(self, *, channel=None, mute=None):
        before = SimpleNamespace(
            channel=self.voice.channel,
            mute=self.voice.mute
        )
        if before.channel and self in before.channel.members:
            before.channel.members.remove(self)
        if channel:
            if self not in channel.members:
                channel.members.append(self)
        self.voice.channel = channel
        if mute is not None:
            self.voice.mute = mute
        after = SimpleNamespace(
            channel=self.voice.channel,
            mute=self.voice.mute
        )
        self.guild.on_voice_state_update(self, before, after)
        return before, after


    return type(
        'MockMember',
        (discord.Member,),
        {
            'bot': bot,
            'display_avatar': SimpleNamespace(url="https://example.com"),
            'display_name': name,
            'edit': edit,
            'guild': guild,
            'id': id,
            'mention': f'<@{id}>',
            'update_voice_state': update_voice_state,
            'voice': SimpleNamespace(channel=channel, mute=False)
        }
    )()

def create_message(allowed_mentions=None, author=None, content=None, channel=None, embeds=None, guild=None, id=None, **kwargs):

    async def add_reaction(self, emoji):
        self.reactions.append(emoji)
    
    async def clear_reactions(self):
        self.reactions.clear()

    async def remove_reaction(self, emoji, user):
        if emoji in self.reactions:
            self.reactions.remove(emoji)
    
    async def edit(self, *, content=None, embed=None, embeds=None, view=None, **kwargs):
        if content is not None:
            self.content = content
        if embed is not None:
            self.edited_embeds.append(embed)
        if embeds is not None:
            self.embeds = embeds
        if view is not None:
            self.view = view
        return self
    
    return type(
        'MockMessage',
        (discord.Message,),
        {
            'add_reaction': add_reaction,
            'allowed_mentions': allowed_mentions,
            'attachments': [],
            'author': author,
            'content': content,
            'channel': channel,
            'clear_reactions': clear_reactions,
            'created_at': datetime.now(timezone.utc),
            'edit': edit,
            'edited_embeds': [],
            'embeds': embeds or [],
            'guild': guild,
            'id': id,
            'reactions': [],
            'remove_reaction': remove_reaction,
            '_state': create_http()
        }
    )()

def create_guild(bot, channels=None, id=None, name=None, members=None, owner_snowflake=None, roles=None):

    def on_voice_state_update(self, member, before, after):
        self.me.dispatch('voice_state_update', member, before, after)

    def get_channel(self, channel_snowflake):
        if channel_snowflake is None:
            return None
        return self.channels.get(channel_snowflake)

    def get_member(self, member_snowflake):
        if member_snowflake is None:
            return None
        return self.members.get(member_snowflake)

    def get_role(self, role_snowflake):
        if role_snowflake is None:
            return None
        return self.roles.get(role_snowflake)

    guild = type(
        'MockGuild',
        (discord.Guild,),
        {
            'id': id,
            'channels': channels,
            'on_voice_state_update': on_voice_state_update,
            'get_channel': get_channel,
            'get_member': get_member,
            "get_role": get_role,
            'me': bot,
            'members': members,
            'name': name,
            'owner_id': owner_snowflake,
            "roles": roles,
            '_state': create_http()
        }
    )()
    return guild


def create_channel(channel_type=None, guild=None, id=None, name=None, object_channel=discord.VoiceChannel):
    
    _messages = []

    def append_message(self, msg):
        _messages.append(msg)

    async def fetch_message(self, message_snowflake):
        if message_snowflake is None:
            return None
        for msg in self.messages:
            if msg.id == message_snowflake:
                return msg
        raise ValueError(f"Message with snowflake {message_snowflake} not found")

    def permissions_for(self, member):
        return SimpleNamespace(send_messages=True)

    async def send_message(author=None, content=None, embeds=None, params=None, **kwargs):
        return {
            'author': {
                'id': author.id,
                'username': author.name
            },
            'channel_id': self.id,
            'content': content,
            'embeds': embeds,
            'id': MESSAGE_ID,
            'params': params
        }

    async def set_permissions(self, target, **overwrites):
        self.overwrites[target.id] = overwrites
        return True

    return type(
        'MockChannel',
        (object_channel,),
        {
            'append_message': append_message,
            'guild': guild,
            'id': id,
            'fetch_message': fetch_message,
            'members': [],
            'mention': f'<#{id}>',
            'messages': _messages,
            'name': name,
            'overwrites': {},
            'permissions_for': permissions_for,
            'set_permissions': set_permissions,
            'send_messages': True,
            'type': channel_type,
            '_state': create_http()
        }
    )()
    
def create_role(guild=None, id=None, name=None, members=None):
    return SimpleNamespace(
        guild=guild,
        id=id,
        name=name,
        members=members,
        mention=f"<@&{id}>"
    )
    
    # channels = {}
    # for cid, (name, channel_type) in channel_defs.items():
    #     channels[cid] = make_mock_channel(channel_type=channel_type, guild=guild, id=cid, name=name)
    # for member in (members or {}).values():
    #     client_channel = list(channels.values())[0]
    #     member.voice.channel = client_channel
    #     client_channel.members.append(member)
    #     for role in roles.values():
    #         role.members.append(member)
    # guild.add_member = lambda member: guild.members.update({member.id: member})


# def make_mock_interaction(bot, guild, channel, user, command_name, options=None):
#     command = bot.tree.get_command(command_name)
#     assert command is not None
#     assert command.id is not None

#     class MockResponse:
#         async def send_message(self, *args, **kwargs): pass
#         async def defer(self, *args, **kwargs): pass
#         def is_done(self): return False

#     class MockFollowup:
#         async def send(self, *args, **kwargs): pass

#     return SimpleNamespace(
#         type=discord.InteractionType.application_command,
#         data={
#             'id': command.id,
#             'name': command_name,
#             'type': 1,
#             'options': [
#                 {'name': k, 'type': 3, 'value': v}
#                 for k, v in (options or {}).items()
#             ]
#         },
#         guild=guild,
#         channel=channel,
#         user=user,
#         client=bot,
#         response=MockResponse(),
#         followup=MockFollowup(),
#         created_at=datetime.now(timezone.utc)
#     )