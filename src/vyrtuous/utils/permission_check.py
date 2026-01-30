import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.guild_dictionary import (
    flush_page,
)
from vyrtuous.inc.helpers import TARGET_PERMISSIONS


class PermissionCheck:

    @classmethod
    async def check_permissions(cls, object_dict, snowflake_kwargs, target):
        bot = DiscordBot.get_instance()
        channel_snowflake = snowflake_kwargs.get("channel_snowflake", None)
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        guild = bot.get_guild(guild_snowflake)
        chunk_size, field_count, lines, pages = 7, 0, [], []
        guild_dictionary = {}
        title = f"{get_random_emoji()} {bot.user.display_name} Missing Permissions"

        if target and target.lower() == "all":
            channel_objs = [
                channel_obj for guild in bot.guilds for channel_obj in guild.channels
            ]
        elif hasattr(object_dict.get("object", None), "channels"):
            channel_objs = object_dict.get("object", None).channels
        else:
            channel_objs = [object_dict.get("object", None)]
        for channel in channel_objs:
            permissions = channel.permissions_for(guild.me)
            missing = []
            for permission in TARGET_PERMISSIONS:
                if not getattr(permissions, permission):
                    missing.append(permission)
            if not missing:
                continue
            guild_dictionary.setdefault(channel.guild.id, {"channels": {}})
            guild_dictionary[channel.guild.id]["channels"].setdefault(channel.id, {})
            guild_dictionary[channel.guild.id]["channels"][channel.id].update(
                {"permissions": missing}
            )
        for guild_snowflake, guild_data in guild_dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                lines.append(f"Channel: {channel.mention}")
                for section_name, permissions in channel_data.items():
                    for permission in permissions:
                        lines.append(f"  ↳ {permission}")
                field_count += 1
                if field_count >= chunk_size:
                    embed.add_field(
                        name="Information", value="\n".join(lines), inline=False
                    )
                    embed, field_count = flush_page(embed, pages, title, guild.name)
                    lines = []
            if lines:
                embed.add_field(
                    name="Information", value="\n".join(lines), inline=False
                )
            pages.append(embed)
        return pages
