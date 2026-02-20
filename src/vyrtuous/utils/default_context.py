from dataclasses import dataclass, field
from typing import Union

import discord
from discord.ext import commands


@dataclass(frozen=True)
class DefaultContext:
    author: discord.Member = field(init=False)
    channel: discord.abc.GuildChannel = field(init=False)
    ctx: commands.Context | None = None
    guild: discord.Guild = field(init=False)
    interaction: discord.Interaction | None = None
    message: discord.Message | None = None
    _source: Union[commands.Context, discord.Interaction, discord.Message, None] = (
        field(init=False)
    )

    def __post_init__(self):
        object.__setattr__(
            self, "_source", self.ctx or self.interaction or self.message
        )
        source = self._source
        if source is None:
            raise ValueError(
                "Must provide at least one of ctx, interaction, or message"
            )
        if isinstance(source, discord.Interaction):
            object.__setattr__(self, "author", source.user)
            object.__setattr__(self, "channel", source.channel)
            object.__setattr__(self, "guild", source.guild)
        else:
            object.__setattr__(self, "author", source.author)
            object.__setattr__(self, "channel", source.channel)
            object.__setattr__(self, "guild", source.guild)

    def to_dict(self) -> dict:
        return {
            "channel_snowflake": self.channel.id,
            "guild_snowflake": self.guild.id,
            "member_snowflake": self.author.id,
        }
