import datetime
from collections import defaultdict

import discord
import pytz
from discord.ext import commands, tasks
from vyrtuous.utils.handlers.mute_service import MuteService
from vyrtuous.utils.handlers.sql_manager import perform_backup, setup_backup_directory
from vyrtuous.utils.inc.helpers import *


class Ruderalis(commands.Cog):

    def __init__(self, bot):
        self.backup_database.start()
        self.bot = bot
        self.config = bot.config
        self.mute_service = MuteService(self.bot.db_pool)

    @commands.after_invoke
    async def after_invoke(self, ctx):
        if hasattr(bot, 'db_pool'):
            await bot.db_pool.close()

    @tasks.loop(hours=24)
    async def backup_database(self):
        try:
            backup_dir = setup_backup_directory('./backups')
            backup_file = perform_backup(
                db_user='postgres',
                db_name='vyrtuous',
                db_host='localhost',
                backup_dir=backup_dir
            )
            logger.info(f'Backup completed successfully: {backup_file}')
        except Exception as e:
            logger.error(f'Error during database backup: {e}')

    @backup_database.before_loop
    async def before_backup(self):
        await self.bot.wait_until_ready()

async def setup(bot: commands.Bot):
    await bot.add_cog(Ruderalis(bot))

