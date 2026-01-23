from discord.ext import commands
import discord


def resolve_author(source):
    if isinstance(source, discord.Interaction):
        member = source.user
    elif isinstance(source, (commands.Context, discord.Message)):
        member = source.author
    return member
