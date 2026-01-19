from discord.ext import commands
import discord


def get_member_snowflake(source):
    if isinstance(source, discord.Interaction):
        member_snowflake = source.user.id
    elif isinstance(source, (commands.Context, discord.Message)):
        member_snowflake = source.author.id
    return member_snowflake
