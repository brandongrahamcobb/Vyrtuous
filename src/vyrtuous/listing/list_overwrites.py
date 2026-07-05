import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.messaging import emojis


def build_embed(obj) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    member_count, role_count, total_count = 0, 0, 0
    obj_name = "All Servers"
    if obj is not None and not isinstance(obj, (int, str)):
        obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Overwrites for {obj_name}"
    if obj == "all":
        for guild in bot.guilds:
            for channel in guild.channels:
                for target_obj, overwrite in channel.overwrites.items():
                    if any(v is not None for v in overwrite):
                        total_count += 1
                        if isinstance(target_obj, discord.Member):
                            member_count += 1
                        else:
                            role_count += 1
    elif isinstance(obj, discord.Guild):
        for channel in obj.channels:
            for target_obj, overwrite in channel.overwrites.items():
                if any(v is not None for v in overwrite):
                    total_count += 1
                    if isinstance(target_obj, discord.Member):
                        member_count += 1
                    else:
                        role_count += 1
    elif isinstance(
        obj, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
    ):
        for target_obj, overwrite in obj.overwrites.items():
            if any(v is not None for v in overwrite):
                total_count += 1
                if isinstance(target_obj, discord.Member):
                    member_count += 1
                else:
                    role_count += 1
    embed = discord.Embed(title=title, color=discord.Color.blue())
    embed.add_field(name="Role overwrites", value=str(role_count), inline=False)
    embed.add_field(name="Member overwrites", value=str(member_count), inline=False)
    embed.add_field(name="Total overwrites", value=str(total_count), inline=False)
    return embed
