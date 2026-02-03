"""text_mute.py The purpose of this program is to inherit from DatabaseFactory to provide the text mute moderation.

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

from vyrtuous.db.base.alias.alias import Alias
from vyrtuous.db.infractions.tmute.text_mute_service import TextMuteService


class TextMuteAlias(Alias):

    category = "tmute"
    service = TextMuteService

    ARGS_MAP = {"alias_name": 1, "member": 2, "duration": 3, "reason": 4}
