"""!/bin/python3
mock_discord_state.py The purpose of this program is to support integration testing for Vyrtuous.

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

import asyncio

from discord.http import HTTPClient
from discord.state import ConnectionState


class MockState(ConnectionState):

    def dispatch(event, *args):
        return None

    handlers = {}
    hooks = {}
    http = HTTPClient(dispatch)
    http._global_over = asyncio.Event()
    http._global_over.set()

    def __init__(self):
        super().__init__(
            dispatch=self.dispatch,
            handlers=self.handlers,
            hooks=self.hooks,
            http=self.http,
        )
        self.channels = {}
