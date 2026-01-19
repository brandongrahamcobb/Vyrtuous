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
        super().__init__(dispatch=self.dispatch, handlers=self.handlers, hooks=self.hooks, http=self.http)
        self._channels = {}