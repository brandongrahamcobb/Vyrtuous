from typing import Union

import discord
from discord.ext import commands

from vyrtuous.commands.author import resolve_author


def source_to_snowflakes(
    source: Union[commands.Context, discord.Interaction, discord.Message],
):
    channel_snowflake = source.channel.id
    guild_snowflake = source.guild.id
    member_snowflake = resolve_author(source=source)
    return {
        "channel_snowflake": channel_snowflake,
        "guild_snowflake": guild_snowflake,
        "member_snowflake": member_snowflake,
    }
