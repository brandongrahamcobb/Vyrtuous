from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot


class SnowflakeContext:
    def __init__(
        self,
        *,
        channel_snowflake: int,
        guild_snowflake: int,
        member_snowflake: int,
    ):
        self.bot: DiscordBot = DiscordBot.get_instance()
        self.guild_snowflake = guild_snowflake
        self.guild = self.bot.get_guild(guild_snowflake)
        if self.guild is None:
            raise GuildNotFound(str(guild_snowflake))
        self.channel_snowflake = channel_snowflake
        self.channel = self.guild.get_channel(channel_snowflake)
        if self.channel is None:
            raise commands.ChannelNotFound(str(channel_snowflake))
        self.member_snowflake = member_snowflake
        self.member = self.guild.get_member(member_snowflake)
