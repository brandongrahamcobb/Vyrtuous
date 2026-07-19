import discord

from vyrtuous.utils.messaging import emojis


def build_embed(obj) -> discord.Embed:
    member_count, role_count, total_count = 0, 0, 0
    obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Overwrites for {obj_name}"
    if isinstance(obj, discord.Guild):
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
