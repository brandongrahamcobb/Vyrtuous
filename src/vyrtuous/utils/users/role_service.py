"""!/bin/python3
role_service.py The purpose of this program is to extend AliasService to service the role class.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import discord
from discord.ext import commands

from vyrtuous.aliases import role_alias_service, unrole_alias_service
from vyrtuous.aliases.alias_context import AliasContext
from vyrtuous.bot.discord_bot import DiscordBot


async def enforce_or_undo(
    alias_ctx: AliasContext, message: discord.Message
) -> discord.Embed:
    bot: DiscordBot = DiscordBot.get_instance()
    guild = bot.get_guild(alias_ctx.guild_snowflake)
    if guild is None:
        raise commands.GuildNotFound(str(alias_ctx.guild_snowflake))
    role = guild.get_role(alias_ctx.role_snowflake)
    if role is None:
        raise commands.RoleNotFound(str(alias_ctx.role_snowflake))
    member = guild.get_member(alias_ctx.member_snowflake)
    if member is None:
        raise commands.MemberNotFound(str(alias_ctx.member_snowflake))
    if role in member.roles:
        unrole_ctx = unrole_alias_service.UnroleMessageContext(
            author_snowflake=message.author.id,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            role_snowflake=alias_ctx.role_snowflake,
        )
        embed = await unrole_alias_service.unrole_by_message(
            ctx=unrole_ctx, display=True
        )
        return embed
    else:
        enrole_ctx = role_alias_service.EnroleMessageContext(
            author_snowflake=message.author.id,
            guild_snowflake=alias_ctx.guild_snowflake,
            member_snowflake=alias_ctx.member_snowflake,
            message_snowflake=message.id,
            message_channel_snowflake=message.channel.id,
            role_snowflake=alias_ctx.role_snowflake,
        )
        embed = await role_alias_service.enrole_by_message(ctx=enrole_ctx, display=True)
        return embed


async def added_role(
    category_class,
    category_role_class,
    guild_snowflake,
    member_snowflake,
    role_snowflake,
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    kwargs = {
        "guild_snowflake": int(guild_snowflake),
        "role_snowflake": str(role_snowflake),
    }
    role = await category_role_class.select(singular=True, **kwargs)
    if role:
        kwargs.update({"channel_snowflake": role.channel_snowflake})
        kwargs.update({"member_snowflake": role.member_snowflake})
        msg = f"Member ({member_snowflake}) was granted the role ({role_snowflake}) for category ({category_class.__name__()}) related to channel ({role.channel_snowflake}) in guild ({guild_snowflake})."
        category = category_class(**kwargs)
        await category.create()
        bot.logger.info(msg)
    else:
        return


async def removed_role(
    category_class,
    category_role_class,
    guild_snowflake,
    member_snowflake,
    role_snowflake,
) -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    kwargs = {
        "guild_snowflake": int(guild_snowflake),
        "role_snowflake": str(role_snowflake),
    }
    role = await category_role_class.select(singular=True, **kwargs)
    if role:
        kwargs.update({"channel_snowflake": role.channel_snowflake})
        kwargs.update({"member_snowflake": role.member_snowflake})
        msg = f"Member ({member_snowflake}) was revoked the role ({role_snowflake}) for category ({category_class.__name__()}) related to channel ({role.channel_snowflake}) in guild ({guild_snowflake})."
        await category_class.delete(**kwargs)
        bot.logger.info(msg)
    else:
        return
