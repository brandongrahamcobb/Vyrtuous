from lucy.utils.setup_logging import logger

import discord
import logging

class Predicator:

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config

    def at_home(self):
        async def predicate(ctx):
            return ctx.guild is not None and ctx.guild.id == self.bot.config['discord_testing_guild_id']
        return commands.check(predicate)

    async def is_vegan(self, user: discord.User):
        async def predicate(ctx):
            guilds = [
                await self.bot.fetch_guild(self.config['discord_testing_guild_id']),
                await self.bot.fetch_guild(730907954345279591)
            ]
            for guild in guilds:
                vegan_role = get(guild.roles, name="Vegan")
                if vegan_role in user.roles:
                    return True
            return False
        return commands.check(predicate)

    def release_mode(self):
        async def predicate(ctx):
            logger.info(f"Release mode setting: {self.bot.config['discord_release_mode']}")
            return ctx.author.id == 154749533429956608 or self.bot.config['discord_release_mode'] or isinstance(ctx.message.channel, discord.DMChannel)
        return commands.check(predicate)

