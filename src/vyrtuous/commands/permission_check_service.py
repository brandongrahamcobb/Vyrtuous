import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import CHUNK_SIZE, TARGET_PERMISSIONS
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_channels,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji


class PermissionCheckService:

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
                PermissionCheckService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
            if skipped_guilds:
                PermissionCheckService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, is_at_home, channel_objs, snowflake_kwargs):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        channel_snowflake = snowflake_kwargs.get("channel_snowflake", None)
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        guild = bot.get_guild(guild_snowflake)
        title = f"{get_random_emoji()} {bot.user.display_name} Missing Permissions"

        dictionary = await PermissionCheckService.build_clean_dictionary(
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
                PermissionCheckService.lines.append(f"Channel: {channel.mention}")
                for section_name, permissions in channel_data.items():
                    for permission in permissions:
                        PermissionCheckService.lines.append(f"  ↳ {permission}")
                field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name="Information",
                        value="\n".join(PermissionCheckService.lines),
                        inline=False,
                    )
                    embed = flush_page(
                        embed, PermissionCheckService.pages, title, guild.name
                    )
                    PermissionCheckService.lines = []
            if PermissionCheckService.lines:
                embed.add_field(
                    name="Information",
                    value="\n".join(PermissionCheckService.lines),
                    inline=False,
                )
            PermissionCheckService.pages.append(embed)
        return PermissionCheckService.pages
