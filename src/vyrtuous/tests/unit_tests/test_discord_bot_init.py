import pytest

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.config import Config
from vyrtuous.db.database import Database


@pytest.mark.asyncio
async def test_discord_bot_init():
    config = Config().get_config()
    db_pool = await Database(config=config).database_init()
    try:
        discord_bot = DiscordBot.get_instance()
    except RuntimeError as e:
        assert isinstance(e, RuntimeError)
    discord_bot = DiscordBot(config=config, db_pool=db_pool)
    await discord_bot.setup_hook()
    # try:
    #     # discord_bot = DiscordBot(config=None, db_pool=db_pool)
    # except Exception as e:
    #     assert isinstance(e, Exception)
    assert discord_bot
