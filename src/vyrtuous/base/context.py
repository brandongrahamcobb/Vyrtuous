from typing import Dict

import discord
from discord.ext import commands

from vyrtuous.commands.author import resolve_author


class Context:
    def __init__(
        self,
        *,
        ctx: commands.Context | None = None,
        interaction: discord.Interaction | None = None,
        message: discord.Message | None = None,
    ):
        self.source_kwargs: Dict[str, int] = {}
        self.source = ctx or interaction or message
        self.author = resolve_author(self.source)
        self.source_channel_snowflake = self.source.channel.id
        self.source_guild_snowflake = self.source.guild.id
        self.source_member_snowflake = self.author.id

    def build_source_kwargs(self):
        self.source_kwargs = {
            "channel_snowflake": self.source_channel_snowflake,
            "guild_snowflake": self.source_guild_snowflake,
            "member_snowflake": self.source_member_snowflake,
        }
