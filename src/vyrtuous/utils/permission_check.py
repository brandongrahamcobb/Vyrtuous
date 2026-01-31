import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.utils.dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_dictionary,
    flush_page,
)
from vyrtuous.inc.helpers import TARGET_PERMISSIONS, CHUNK_SIZE


class PermissionCheck:

    lines, pages = [], []

    @classmethod
    async def build_clean_dictionary(cls, channel_objs, is_at_home, me):
        dictionary = {}
        for channel in channel_objs:
            permissions = channel.permissions_for(me)
            missing = []
            for permission in TARGET_PERMISSIONS:
                if not getattr(permissions, permission):
                    missing.append(permission)
            if not missing:
                continue
            dictionary.setdefault(channel.guild.id, {"channels": {}})
            dictionary[channel.guild.id]["channels"].setdefault(channel.id, {})
            dictionary[channel.guild.id]["channels"][channel.id].update(
                {"permissions": missing}
            )
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_channels:
                PermissionCheck.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
            if skipped_guilds:
                PermissionCheck.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, is_at_home, channel_objs, snowflake_kwargs):
        bot = DiscordBot.get_instance()
        channel_snowflake = snowflake_kwargs.get("channel_snowflake", None)
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        guild = bot.get_guild(guild_snowflake)
        title = f"{get_random_emoji()} {bot.user.display_name} Missing Permissions"

        dictionary = await PermissionCheck.build_clean_dictionary(
            channel_objs=channel_objs, is_at_home=is_at_home, me=guild.me
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, channel_data in guild_data.get(
                "channels", {}
            ).items():
                channel = guild.get_channel(channel_snowflake)
                PermissionCheck.lines.append(f"Channel: {channel.mention}")
                for section_name, permissions in channel_data.items():
                    for permission in permissions:
                        PermissionCheck.lines.append(f"  ↳ {permission}")
                field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(PermissionCheck.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, PermissionCheck.pages, title, guild.name)
                    PermissionCheck.lines = []
            if PermissionCheck.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(PermissionCheck.lines),
                    inline=False,
                )
            PermissionCheck.pages.append(embed)
        return PermissionCheck.pages
