"""room.py The purpose of this program is to inherit from the DatabaseFactory to provide a parent to room classes.

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

from vyrtuous.database.actions.alias import Alias
from vyrtuous.database.actions.ban import Ban
from vyrtuous.database.actions.flag import Flag
from vyrtuous.database.actions.text_mute import TextMute
from vyrtuous.database.actions.voice_mute import VoiceMute
from vyrtuous.database.logs.history import History
from vyrtuous.database.roles.coordinator import Coordinator
from vyrtuous.database.roles.moderator import Moderator
from vyrtuous.database.roles.vegan import Vegan
from vyrtuous.database.rooms.stage import Stage
from vyrtuous.database.rooms.temporary_room import TemporaryRoom
from vyrtuous.database.rooms.video_room import VideoRoom
from vyrtuous.database.database_factory import DatabaseFactory

member_relevant_objects_dict = {
    Ban.SINGULAR: Ban,
    Coordinator.SINGULAR: Coordinator,
    Flag.SINGULAR: Flag,
    Moderator.SINGULAR: Moderator,
    TextMute.SINGULAR: TextMute,
    Vegan.SINGULAR: Vegan,
    VoiceMute.SINGULAR: VoiceMute,
}

room_relevant_objects_dict = {
    Alias.SINGULAR: Alias,
    Ban.SINGULAR: Ban,
    Coordinator.SINGULAR: Coordinator,
    Flag.SINGULAR: Flag,
    History.SINGULAR: History,
    Moderator.SINGULAR: Moderator,
    TemporaryRoom.SINGULAR: TemporaryRoom,
    TextMute.SINGULAR: TextMute,
    Vegan.SINGULAR: Vegan,
    VideoRoom.SINGULAR: VideoRoom,
    VoiceMute.SINGULAR: VoiceMute,
    Stage.SINGULAR: Stage,
}


class Room(DatabaseFactory):
    pass
