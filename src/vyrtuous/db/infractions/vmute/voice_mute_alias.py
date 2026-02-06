"""!/bin/python3
voice_mute_alias.py The purpose of this program is to extend Alias to provide the voice mute infraction alias.

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

from vyrtuous.db.alias.alias import Alias
from vyrtuous.db.infractions.vmute.voice_mute_service import VoiceMuteService


class VoiceMuteAlias(Alias):

    category = "vmute"
    service = VoiceMuteService

    ARGS_MAP = {"alias_name": 1, "member": 2, "duration": 3, "reason": 4}
