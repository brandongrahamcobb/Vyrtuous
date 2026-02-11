"""!/bin/python3
flag_service.py The purpose of this program is to extend AliasService to service flag infractions.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord

from vyrtuous.base.record_service import RecordService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.flag.flag import Flag
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.stream.stream_service import StreamService
from vyrtuous.utils.dictionary import (
    clean_dictionary,
    flush_page,
    generate_skipped_dict_pages,
    generate_skipped_guilds,
    generate_skipped_members,
    generate_skipped_set_pages,
)
from vyrtuous.utils.emojis import get_random_emoji


class FlagService(RecordService):
    lines, pages = [], []
    model = Flag

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        dictionary = {}
        flags = await Flag.select(singular=False, **where_kwargs)
        for flag in flags:
            dictionary.setdefault(flag.guild_snowflake, {"members": {}})
            dictionary[flag.guild_snowflake]["members"].setdefault(
                flag.member_snowflake, {"flags": {}}
            )
            dictionary[flag.guild_snowflake]["members"][flag.member_snowflake][
                "flags"
            ].setdefault(flag.channel_snowflake, {})
            dictionary[flag.guild_snowflake]["members"][flag.member_snowflake]["flags"][
                flag.channel_snowflake
            ].update(
                {
                    "reason": flag.reason,
                }
            )
        skipped_guilds = generate_skipped_guilds(dictionary)
        skipped_members = generate_skipped_members(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_guilds=skipped_guilds,
            skipped_members=skipped_members,
        )
        if is_at_home:
            if skipped_guilds:
                FlagService.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_members:
                FlagService.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_members,
                        title="Skipped Members in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        cls.lines = []
        cls.pages = []
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Flags {f'for {object_dict.get('name', None)}' if isinstance(object_dict.get('object', None), discord.Member) else ''}"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await FlagService.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        flag_n = 0
        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            thumbnail = False
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for member_snowflake, flag_dictionary in guild_data.get("members").items():
                member = guild.get_member(member_snowflake)
                if not member:
                    continue
                if not isinstance(object_dict.get("object", None), discord.Member):
                    FlagService.lines.append(
                        f"**User:** {member.display_name} {member.mention}"
                    )
                elif not thumbnail:
                    embed.set_thumbnail(
                        url=object_dict.get("object", None).display_avatar.url
                    )
                    thumbnail = True
                for channel_snowflake, channel_dictionary in flag_dictionary.get(
                    "flags", {}
                ).items():
                    channel = guild.get_channel(channel_snowflake)
                    if not isinstance(
                        object_dict.get("object"), discord.abc.GuildChannel
                    ):
                        FlagService.lines.append(f"**Channel:** {channel.mention}")
                    if isinstance(object_dict.get("object"), discord.Member):
                        FlagService.lines.append(
                            f"**Reason:** {channel_dictionary['reason']}"
                        )
                    flag_n += 1
                    field_count += 1
                    if field_count >= CHUNK_SIZE:
                        embed.add_field(
                            name="Information",
                            value="\n".join(FlagService.lines),
                            inline=False,
                        )
                        embed = flush_page(embed, FlagService.pages, title, guild.name)
                        FlagService.lines = []
            if FlagService.lines:
                embed.add_field(
                    name="Information", value="\n".join(FlagService.lines), inline=False
                )
            FlagService.pages.append(embed)
        if FlagService.pages:
            FlagService.pages[0].description = f"**({flag_n})**"
        return FlagService.pages

    @classmethod
    async def enforce(cls, ctx, source, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        flag = Flag(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
            reason=ctx.reason,
        )
        await flag.create()
        bot = DiscordBot.get_instance()
        cog = bot.get_cog("ChannelEventListeners")
        cog.flags.append(flag)
        await StreamService.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="flag",
            member=member,
            source=source,
            reason=ctx.reason,
        )
        embed = await FlagService.act_embed(ctx=ctx)
        return await state.end(success=embed)

    @classmethod
    async def undo(cls, ctx, source, state):
        bot = DiscordBot.get_instance()
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        await Flag.delete(
            channel_snowflake=ctx.target_channel_snowflake,
            guild_snowflake=ctx.source_guild_snowflake,
            member_snowflake=ctx.target_member_snowflake,
        )
        bot = DiscordBot.get_instance()
        cog = bot.get_cog("ChannelEventListeners")
        for flag in cog.flags:
            if flag.channel_snowflake == ctx.target_channel_snowflake:
                cog.flags.remove(flag)
                break
        await StreamService.send_entry(
            channel_snowflake=ctx.target_channel_snowflake,
            identifier="unflag",
            is_modification=True,
            member=member,
            source=source,
        )
        embed = await FlagService.undo_embed(ctx=ctx)
        return await state.end(success=embed)

    @classmethod
    async def act_embed(cls, ctx):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(ctx.target_channel_snowflake)
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"{get_random_emoji()} {member.display_name} has been flagged",
            description=(
                f"**User:** {member.mention}\n"
                f"**Channel:** {channel.mention}\n"
                f"**Reason:** {ctx.reason}"
            ),
            color=discord.Color.blue(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed

    @classmethod
    async def undo_embed(cls, ctx):
        bot = DiscordBot.get_instance()
        channel = bot.get_channel(ctx.target_channel_snowflake)
        guild = bot.get_guild(ctx.source_guild_snowflake)
        member = guild.get_member(ctx.target_member_snowflake)
        embed = discord.Embed(
            title=f"{get_random_emoji()} {member.display_name} has been unflagged",
            description=(f"**User:** {member.mention}\n**Channel:** {channel.mention}"),
            color=discord.Color.yellow(),
        )
        embed.set_thumbnail(url=member.display_avatar.url)
        return embed
