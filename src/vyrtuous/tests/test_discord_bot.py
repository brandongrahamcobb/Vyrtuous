import pytest
from discord.ext import commands

@pytest.mark.asyncio
async def make_mock_context(bot: commands.Bot, command_name: str):
    channel = type("Channel", (), {"id": 1234567890, "send": lambda *args, **kwargs: None})()
    author = type("User", (), {"id": 9999, "name": "TestUser"})()
    message = type("Message", (), {"content": f"!{command_name}", "author": author, "channel": channel})()
    ctx = await bot.get_context(message)
    return ctx