# """make_mock_objects.py The purpose of this program is to provide the shared test variables for tests using Discord objects.
# Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# """

# from datetime import datetime, timezone
# from types import SimpleNamespace

# from discord.http import HTTPClient
# from discord.state import ConnectionState
# import discord

# from vyrtuous.inc.helpers import NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE, NOT_PRIVILEGED_AUTHOR_NAME_ONE


# def create_state():
#     dispatch = lambda event, *args: None
#     handlers = {}
#     hooks = {}
#     http = HTTPClient(dispatch)
#     # http._global_over = asyncio.Event()
#     # http._global_over.set()
#     state = ConnectionState(
#         dispatch=dispatch, handlers=handlers, hooks=hooks, http=http
#     )
#     state._channels = {}
#     return state


# def create_member(bot=False, channel=None, guild=None, id=None, name=None):

#     data = {
#         "user": {
#             "id": id,  # member's user ID
#             "username": NOT_PRIVILEGED_AUTHOR_SNOWFLAKE_ONE,  # username
#             "discriminator": "1234",  # discriminator
#             "avatar": None,  # avatar hash
#             "bot": False,  # whether the user is a bot
#         },
#         "nick": NOT_PRIVILEGED_AUTHOR_NAME_ONE,  # nickname in the guild
#         "roles": [],  # list of role IDs
#         "joined_at": "2025-01-01T00:00:00.000000+00:00",  # ISO timestamp
#         "premium_since": None,  # boosting start date
#         "pending": False,  # if membership is pending verification
#         "deaf": False,  # voice state
#         "mute": False,  # voice state
#         "communication_disabled_until": None,  # timeout until
#         "flags": 0,  # member flags
#         "avatar": None,
#     }

#     async def edit(self, **kwargs):
#         for k, v in kwargs.items():
#             setattr(self, k, v)
#         return self

#     async def update_voice_state(self, *, channel=None, mute=None):
#         before = SimpleNamespace(channel=self.voice.channel, mute=self.voice.mute)
#         if before.channel and self in before.channel.members:
#             before.channel.members.remove(self)
#         if channel:
#             if self not in channel.members:
#                 channel.members.append(self)
#         self.voice.channel = channel
#         if mute is not None:
#             self.voice.mute = mute
#         after = SimpleNamespace(channel=self.voice.channel, mute=self.voice.mute)
#         self.guild.on_voice_state_update(self, before, after)
#         return before, after

#     return type(
#         "MockMember",
#         (discord.Member,),
#         {
#             "bot": bot,
#             "display_avatar": SimpleNamespace(url="https://example.com"),
#             "display_name": name,
#             "edit": edit,
#             "guild": guild,
#             "id": id,
#             "mention": f"<@{id}>",
#             "update_voice_state": update_voice_state,
#             "voice": SimpleNamespace(channel=channel, mute=False),
#         },
#     )(state=create_state(), guild=guild, data=data)


# def create_message(
#     allowed_mentions=None,
#     author=None,
#     bot=False,
#     content=None,
#     channel=None,
#     embed=None,
#     embeds=None,
#     file=None,
#     guild=None,
#     id=None,
#     paginated=None,
#     **kwargs,
# ):

#     data = {
#         "id": id,  # Message ID
#         "channel_id": channel.id,  # The channel ID the message belongs to
#         "guild_id": guild.id,  # Guild ID (optional for DMs)
#         "author": {  # User who sent the message
#             "id": author.id,
#             "username": author.name,
#             "discriminator": "1234",
#             "avatar": None,
#             "bot": False,
#             "system": False,
#             "public_flags": 0,
#         },
#         "member": {  # Optional guild-specific member info
#             "nick": "Tester",
#             "roles": [],
#             "joined_at": "2025-01-01T00:00:00.000000+00:00",
#             "premium_since": None,
#             "deaf": False,
#             "mute": False,
#             "pending": False,
#             "permissions": 104189505,  # example permissions integer
#             "communication_disabled_until": None,
#         },
#         "content": content,  # The message content
#         "timestamp": "2025-01-04T12:00:00.000000+00:00",  # creation timestamp ISO format
#         "edited_timestamp": None,
#         "tts": False,
#         "mention_everyone": False,
#         "mentions": [],  # list of mentioned user objects
#         "mention_roles": [],  # list of mentioned role IDs
#         "attachments": [],  # list of attachments
#         "embeds": embeds if embeds else [],  # list of embeds
#         "reactions": [],  # list of reactions
#         "pinned": False,
#         "webhook_id": None,
#         "type": 0,  # message type, usually 0 (default)
#         "activity": None,
#         "application": {
#             "id": 0,
#             "flags": 0,
#             "cover_image": None,
#             "description": "Test application",
#             "icon": None,
#             "name": "TestApp",
#             "team": None,
#             "verify_key": None,
#         },
#         "message_reference": {
#             "message_id": 123,
#             "channel_id": channel.id,
#             "guild_id": guild.id if guild else None,
#             "type": 0,
#         },
#         "flags": 0,
#     }

#     return type(
#         "MockMessage",
#         (discord.Message,),
#         {
#             "__init__": lambda self: discord.Message.__init__(
#                 self, state=create_state(), channel=channel, data=data
#             ),
#             "add_reaction": add_reaction,
#             "allowed_mentions": allowed_mentions,
#             "attachments": [],
#             "author": author,
#             "content": content,
#             "channel": channel,
#             "clear_reactions": clear_reactions,
#             "created_at": datetime.now(timezone.utc),
#             "edit": edit,
#             "edited_embeds": [],
#             "embed": {},
#             "embeds": [],
#             "file": file,
#             "flags": 0,
#             "guild": guild,
#             "id": id,
#             "reactions": [],
#             "remove_reaction": remove_reaction,
#             "type": 0,
#         },
#     )()


# def create_guild(
#     bot,
#     channels=None,
#     id=None,
#     name=None,
#     members=None,
#     owner_snowflake=None,
#     roles=None,
# ):

#     data = {
#         "id": id,
#         "name": name,
#         "icon": None,
#         "splash": None,
#         "discovery_splash": None,
#         "owner_id": owner_snowflake,
#         "region": None,
#         "afk_channel_id": None,
#         "afk_timeout": 300,
#         "widget_enabled": False,
#         "widget_channel_id": None,
#         "verification_level": 0,
#         "default_message_notifications": 0,
#         "explicit_content_filter": 0,
#         "roles": [],
#         "emojis": [],
#         "features": [],
#         "mfa_level": 0,
#         "application_id": None,
#         "system_channel_id": None,
#         "system_channel_flags": 0,
#         "rules_channel_id": None,
#         "joined_at": "2025-01-01T00:00:00.000000+00:00",
#         "large": False,
#         "member_count": 0,
#         "voice_states": {},
#         "members": [],
#         "channels": [],
#         "presences": {},
#         "max_presences": None,
#         "max_members": None,
#         "vanity_url_code": None,
#         "description": None,
#         "banner": None,
#         "premium_tier": 0,
#         "premium_subscription_count": 0,
#         "preferred_locale": "en-US",
#         "public_updates_channel_id": None,
#         "max_video_channel_users": 25,
#         "approximate_member_count": None,
#         "approximate_presence_count": None,
#         "nsfw_level": 0,
#     }

#     def on_voice_state_update(self, member, before, after):
#         self.me.dispatch("voice_state_update", member, before, after)

#     def get_channel(self, channel_snowflake):
#         if channel_snowflake is None:
#             return None
#         return self.channels.get(channel_snowflake, None)

#     def get_member(self, member_snowflake):
#         if member_snowflake is None:
#             return None
#         for member in self.members:
#             if member.id == member_snowflake:
#                 return member

#     def get_role(self, role_snowflake):
#         if role_snowflake is None:
#             return None
#         for role in self._roles:
#             if role.id == role_snowflake:
#                 return role

#     return type(
#         "MockGuild",
#         (discord.Guild,),
#         {
#             "id": id,
#             "channels": channels,
#             "on_voice_state_update": on_voice_state_update,
#             "get_channel": get_channel,
#             "get_member": get_member,
#             "get_role": get_role,
#             "icon": SimpleNamespace(url="https://example.com"),
#             "me": bot,
#             "members": members,
#             "name": name,
#             "owner_id": owner_snowflake,
#             "roles": roles,
#         },
#     )(data=data, state=create_state())


# def create_channel(
#     channel_type=None,
#     guild=None,
#     id=None,
#     name=None,
#     object_channel=discord.VoiceChannel,
# ):

#     _messages = []
#     data = {
#         "id": id,
#         "name": name,
#         "type": channel_type or 0,
#         "position": 0,  # required
#         "nsfw": False,  # usually required for text channels
#         "topic": "",  # for TextChannel
#         "bitrate": 64000,  # for VoiceChannel
#         "user_limit": 0,  # for VoiceChannel
#         "rate_limit_per_user": 0,
#         "parent_id": None,
#         "guild_id": getattr(guild, "id", None),
#     }

#     def append_message(self, msg):
#         _messages.append(msg)

#     async def fetch_message(self, message_snowflake):
#         if message_snowflake is None:
#             return None
#         for msg in self.messages:
#             if msg.id == message_snowflake:
#                 return msg
#         raise ValueError(f"Message with snowflake {message_snowflake} not found")

#     def permissions_for(self, member):
#         return SimpleNamespace(send_messages=True)

#     async def set_permissions(self, target, **overwrites):
#         self.overwrites[target.id] = overwrites
#         return True

#     return type(
#         "MockChannel",
#         (object_channel,),
#         {
#             "append_message": append_message,
#             "guild": guild,
#             "id": id,
#             "fetch_message": fetch_message,
#             "members": [],
#             "mention": f"<#{id}>",
#             "messages": _messages,
#             "name": name,
#             "overwrites": {},
#             "permissions_for": permissions_for,
#             "set_permissions": set_permissions,
#             # 'send': send,
#             "send_messages": True,
#             "type": channel_type,
#         },
#     )(state=create_state(), guild=guild, data=data)


# def create_role(guild=None, id=None, name=None, members=None, object_role=discord.Role):
#     data = {
#         "id": id,
#         "name": name,
#         "color": 0,
#         "hoist": False,
#         "position": 0,
#         "permissions": 0,
#         "managed": False,
#         "mentionable": False,
#         "guild_id": getattr(guild, "id", None),
#     }

#     def mention(self):
#         return f"<@&{self.id}"

#     return type(
#         "MockRole",
#         (object_role,),
#         {
#             "guild": guild,
#             "id": id,
#             "name": name,
#             "members": members or [],
#             "mention": f"<@&{id}>",
#             "is_default": False,
#             "color": 0,
#             "position": 0,
#             "permissions": 0,
#             "managed": False,
#             "mentionable": False,
#         },
#     )(state=guild._state, guild=guild, data=data)

#     # channels = {}
#     # for cid, (name, channel_type) in channel_defs.items():
#     #     channels[cid] = make_mock_channel(channel_type=channel_type, guild=guild, id=cid, name=name)
#     # for member in (members or {}).values():
#     #     client_channel = list(channels.values())[0]
#     #     member.voice.channel = client_channel
#     #     client_channel.members.append(member)
#     #     for role in roles.values():
#     #         role.members.append(member)
#     # guild.add_member = lambda member: guild.members.update({member.id: member})


# # def make_mock_interaction(bot, guild, channel, user, command_name, options=None):
# #     command = bot.tree.get_command(command_name)
# #     assert command is not None
# #     assert command.id is not None

# #     class MockResponse:
# #         async def send_message(self, *args, **kwargs): pass
# #         async def defer(self, *args, **kwargs): pass
# #         def is_done(self): return False

# #     class MockFollowup:
# #         async def send(self, *args, **kwargs): pass

# #     return SimpleNamespace(
# #         type=discord.InteractionType.application_command,
# #         data={
# #             'id': command.id,
# #             'name': command_name,
# #             'type': 1,
# #             'options': [
# #                 {'name': k, 'type': 3, 'value': v}
# #                 for k, v in (options or {}).items()
# #             ]
# #         },
# #         guild=guild,
# #         channel=channel,
# #         user=user,
# #         client=bot,
# #         response=MockResponse(),
# #         followup=MockFollowup(),
# #         created_at=datetime.now(timezone.utc)
# #     )
