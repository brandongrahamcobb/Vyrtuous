"""!/bin/python3
alias_service.py The purpose of this program is to extend Service to service aliases.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from typing import Union

from vyrtuous.aliases.ban_alias import BanAlias
from vyrtuous.aliases.flag_alias import FlagAlias
from vyrtuous.aliases.role_alias import RoleAlias
from vyrtuous.aliases.text_mute_alias import TextMuteAlias
from vyrtuous.aliases.voice_mute_alias import VoiceMuteAlias
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.alias import Alias, NotAlias
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.utils.errors.error import ChannelNotFound, GuildNotFound, RoleNotFound

MODEL = Alias
ALIASES = [BanAlias, FlagAlias, RoleAlias, TextMuteAlias, VoiceMuteAlias]
CATEGORY_TO_HELP = {
    "ban": [
        "**member**: Specify a member ID/mention.",
        "**duration**: Specify a duration in m/h/d.",
        "**reason**: Specify a reason for the ban.",
    ],
    "flag": [
        "**member**: Specify a member ID/mention.",
        "**reason**: Specify a reason for the flag.",
    ],
    "role": [
        "**member**: Specify a member ID/mention.",
    ],
    "tmute": [
        "**member**: Specify a member ID/mention.",
        "**duration**: Specify a duration in m/h/d.",
        "**reason**: Specify a reason for the text-mute.",
    ],
    "vmute": [
        "**member**: Specify a member ID/mention.",
        "**duration**: Specify a duration in m/h/d.",
        "**reason**: Specify a reason for the voice-mute.",
    ],
}
CATEGORY_TO_DESCRIPTION = {
    "ban": "Toggles a ban.",
    "flag": "Toggles a moderation flag.",
    "role": "Toggles a role to a user.",
    "tmute": "Toggles a mute in text channels.",
    "vmute": "Toggles a mute in voice channels.",
}
CATEGORY_TO_PERMISSION = {
    "ban": "command.moderation.ban",
    "flag": "command.moderation.flag",
    "tmute": "commands.moderation.text-mute",
    "vmute": "commands.moderation.voice-mute",
}


async def disable(alias_name: str, guild_snowflake: int) -> str:
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    alias = await database_factory.select(
        alias_name=alias_name, guild_snowflake=guild_snowflake, singular=True
    )
    if not alias:
        return f"No aliases found for `{alias_name}`."
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(alias.channel_snowflake)
    if channel is None:
        raise ChannelNotFound(alias.channel_snowflake)
    if getattr(alias, "role_snowflake"):
        role = guild.get_role(alias.role_snowflake)
        if not role:
            raise RoleNotFound(alias.role_snowflake)
        msg = (
            f"Alias `{alias.alias_name}` of type "
            f"`{alias.category}` for channel {channel.mention if channel else alias.channel_snowflake} "
            f" and role {role.mention} deleted successfully."
        )
    else:
        msg = (
            f"Alias `{alias.alias_name}` of type "
            f"`{alias.category}` for channel {channel.mention} "
            f"deleted successfully."
        )
    return msg


async def enable(
    alias_name: str,
    category: str,
    channel_snowflake: int,
    guild_snowflake: int,
    role_snowflake: int | None = None,
) -> str:
    if category == "role" and role_snowflake is None:
        return f"Alias of type `{category}` requires a role snowflake."
    bot: DiscordBot = DiscordBot.get_instance()
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    guild = bot.get_guild(guild_snowflake)
    if guild is None:
        raise GuildNotFound(str(guild_snowflake))
    channel = guild.get_channel(channel_snowflake)
    if channel is None:
        raise ChannelNotFound(str(channel_snowflake))
    msg = (
        f"Alias `{alias_name}` of type `{category}` "
        f"created successfully for channel {channel.mention}."
    )
    alias = await database_factory.select(
        category=category,
        alias_name=alias_name,
        guild_snowflake=guild_snowflake,
        singular=True,
    )
    if alias and alias.category != "role":
        return (
            f"Alias of type `{category}` "
            f"with the name `{alias_name} already exists in this server channel ({guild.name})."
        )
    if role_snowflake:
        role = guild.get_role(role_snowflake)
        if role is None:
            raise RoleNotFound(str(role_snowflake))
        msg = (
            f"Alias `{alias_name}` of type `{category}` "
            f"created successfully for channel {channel.mention} with role {role.mention}."
        )
    return msg


def alias_category_to_alias(
    category: str,
) -> Union[BanAlias, FlagAlias, RoleAlias, TextMuteAlias, VoiceMuteAlias]:
    for a in ALIASES:
        if a.category == category:
            return a()
    raise NotAlias()
