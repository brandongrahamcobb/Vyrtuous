"""action.py The purpose of this program is to be a child of DatabaseFactory and the parent to all moderation actions.

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

from datetime import datetime, timezone


import discord

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.fields.duration import DurationObject
from vyrtuous.utils.highest_role import resolve_highest_role
from vyrtuous.utils.dictionary import (
    generate_skipped_dict_pages,
    generate_skipped_set_pages,
    generate_skipped_channels,
    generate_skipped_guilds,
    clean_dictionary,
    flush_page,
)
from vyrtuous.utils.emojis import get_random_emoji
from vyrtuous.inc.helpers import CHUNK_SIZE
from vyrtuous.cogs.aliases import Aliases


class Alias(DatabaseFactory):

    identifier: str

    ACT = "alias"
    PLURAL = "Aliases"
    SINGULAR = "Alias"
    SCOPES = ["channels"]
    UNDO = "xalias"

    REQUIRED_INSTANTIATION_ARGS = [
        "alias_name",
        "category",
        "channel_snowflake",
        "guild_snowflake",
    ]
    OPTIONAL_ARGS = [
        "created_at",
        "expired",
        "expires_in",
        "reason",
        "role_snowflake",
        "updated_at",
    ]

    TABLE_NAME = "command_aliases"
    lines, pages = [], []

    def __init__(
        self,
        alias_name: str,
        category: str,
        channel_snowflake: int,
        guild_snowflake: int,
        created_at: datetime = datetime.now(timezone.utc),
        role_snowflake: int = 0,
        updated_at: datetime = datetime.now(timezone.utc),
        **kwargs,
    ):
        super().__init__()
        self.bot = DiscordBot.get_instance()
        self.alias_name = alias_name
        self.category = category
        self.channel_snowflake = channel_snowflake
        self.channel_mention = f"<#{channel_snowflake}>"
        self.created_at = created_at
        self.guild_snowflake = guild_snowflake
        self.role_snowflake = role_snowflake
        self.role_mention = f"<@&{role_snowflake}>"
        self.updated_at = updated_at

    def for_category(self) -> type['Alias'] | None:
        for subclass in self.__class__.__subclasses__():
            if getattr(subclass, "category", None) == self.category.lower():
                return subclass
        return None
    
    async def build_existing_information(self, infraction_information, member_snowflake):
        kwargs = {}
        cls = self.for_category()
        if cls is None:
            raise ValueError(f"No Alias subclass found for category {self.category!r}.")
        primary_keys = await cls.primary_keys()
        if "channel_snowflake" in primary_keys:
            kwargs.update({"channel_snowflake": self.channel_snowflake})
        if "guild_snowflake" in primary_keys:
            kwargs.update({"guild_snowflake": self.guild_snowflake})
        if "member_snowflake" in primary_keys:
            kwargs.update({"member_snowflake": member_snowflake})
        infraction_existing = await cls.select(
            **kwargs,
            singular=True,
        )
        infraction_modification = False
        if infraction_existing:
            infraction_modification = True
        infraction_information.update(
            {
                "infraction_existing": infraction_existing,
                "infraction_modification": infraction_modification,
                "infraction_role_snowflake": (
                    self.role_snowflake if self.role_snowflake else None
                ),
            }
        )
        return infraction_information

    async def build_duration_information(
        self, infraction_information, duration, state
    ):
        channel = self.bot.get_channel(self.channel_snowflake)
        cap = await Cap.select(
            category=self.category,
            channel_snowflake=str(self.channel_snowflake),
            guild_snowflake=str(self.guild_snowflake),
            singular=True,
        )
        if not hasattr(cap, "duration"):
            cap_duration = DurationObject("8h").to_seconds()
        else:
            cap_duration = cap.duration_seconds
        if (
            duration.to_timedelta().total_seconds() > cap_duration
            or duration.number == 0
        ):
            if infraction_information["executor_role"] == "Moderator":
                duration_str = DurationObject.from_seconds(cap_duration)
                return await state.end(
                    warning=f"Cannot set the "
                    f"{infraction_information['alias_class'].SINGULAR} beyond {duration_str} as a "
                    f"{infraction_information['executor_role']} in {channel.mention}."
                )
        expires_in = datetime.now(timezone.utc) + duration.to_timedelta()
        infraction_information.update(
            {
                "infraction_cap_duration": cap_duration,
                "infraction_duration": duration,
                "infraction_expires_in": expires_in,
            }
        )
        return infraction_information

    async def build_infraction_information(
        self, author_snowflake, duration, infraction_class, member_snowflake, reason, state
    ):
        infraction_information = {}
        if getattr(self, "role_snowflake"):
            infraction_information.update(
                {"infraction_role_snowflake": self.role_snowflake}
            )
        infraction_information.update(
            {
                "infraction_channel_snowflake": self.channel_snowflake,
                "infraction_guild_snowflake": self.guild_snowflake,
                "infraction_member_snowflake": member_snowflake,
            }
        )
        infraction_executor_role = await resolve_highest_role(
            channel_snowflake=self.channel_snowflake,
            guild_snowflake=self.guild_snowflake,
            member_snowflake=author_snowflake,
        )
        infraction_information.update({"infraction_executor_role": infraction_executor_role})
        infraction_information = await self.build_duration_information(
            infraction_information=infraction_information,
            duration=duration,
            state=state,
        )
        infraction_information.update({"infraction_reason": reason})
        infraction_information = await self.build_existing_information(
            infraction_information=infraction_information, member_snowflake=int(member_snowflake)
        )
        return infraction_information

    @classmethod
    async def build_clean_dictionary(cls, is_at_home, where_kwargs):
        bot = DiscordBot.get_instance()
        dictionary = {}
        aliases = await Alias.select(singular=False, **where_kwargs)
        for alias in aliases:
            dictionary.setdefault(alias.guild_snowflake, {"channels": {}})
            dictionary[alias.guild_snowflake]["channels"].setdefault(
                alias.channel_snowflake, {"aliases": {}}
            )
            dictionary[alias.guild_snowflake]["channels"][alias.channel_snowflake][
                "aliases"
            ].setdefault(alias.category, {})[alias.alias_name] = []
            if alias.category == "role":
                guild = bot.get_guild(alias.guild_snowflake)
                if guild:
                    role = guild.get_role(alias.role_snowflake)
                    dictionary[alias.guild_snowflake]["channels"][
                        alias.channel_snowflake
                    ]["aliases"][alias.category][alias.alias_name] = role.mention
        skipped_channels = generate_skipped_channels(dictionary)
        skipped_guilds = generate_skipped_guilds(dictionary)
        cleaned_dictionary = clean_dictionary(
            dictionary=dictionary,
            skipped_channels=skipped_channels,
            skipped_guilds=skipped_guilds,
        )
        if is_at_home:
            if skipped_guilds:
                Alias.pages.extend(
                    generate_skipped_set_pages(
                        skipped=skipped_guilds,
                        title="Skipped Servers",
                    )
                )
            if skipped_channels:
                Alias.pages.extend(
                    generate_skipped_dict_pages(
                        skipped=skipped_channels,
                        title="Skipped Channels in Server",
                    )
                )
        return cleaned_dictionary

    @classmethod
    async def build_pages(cls, object_dict, is_at_home):
        bot = DiscordBot.get_instance()
        title = f"{get_random_emoji()} Command Aliases"

        where_kwargs = object_dict.get("columns", None)
        dictionary = await Alias.build_clean_dictionary(
            is_at_home=is_at_home, where_kwargs=where_kwargs
        )

        for guild_snowflake, guild_data in dictionary.items():
            field_count = 0
            guild = bot.get_guild(guild_snowflake)
            embed = discord.Embed(
                title=title, description=guild.name, color=discord.Color.blue()
            )
            for channel_snowflake, dictionary in guild_data.get("channels", {}).items():
                channel = guild.get_channel(channel_snowflake)
                for category, alias_data in dictionary["aliases"].items():
                    Alias.lines.append(f"{category}")
                    for name, role_mention in alias_data.items():
                        if category == "role":
                            Alias.lines.append(f"  ↳ `{name}` -> {role_mention}")
                        else:
                            Alias.lines.append(f"  ↳ `{name}`")
                    field_count += 1
                if field_count >= CHUNK_SIZE:
                    embed.add_field(
                        name=f"Channel: {channel.mention}",
                        value="\n".join(Alias.lines),
                        inline=False,
                    )
                    embed = flush_page(embed, Alias.pages, title, guild.name)
                    Alias.lines = []
                    field_count = 0
            if Alias.lines:
                embed.add_field(
                    name=f"Channel: {channel.mention}",
                    value="\n".join(Alias.lines),
                    inline=False,
                )
            Alias.pages.append(embed)
        return Alias.pages

    @classmethod
    async def delete_alias(cls, alias_name, snowflake_kwargs):
        bot = DiscordBot.get_instance()
        guild_snowflake = snowflake_kwargs.get("guild_snowflake", None)
        where_kwargs = {
            "alias_name": alias_name,
            "guild_snowflake": int(guild_snowflake),
        }
        alias = await Alias.select(singular=True, **where_kwargs)
        if not alias:
            return f"No aliases found for `{alias_name}`."
        guild = bot.get_guild(guild_snowflake)
        channel = guild.get_channel(alias.channel_snowflake)
        if getattr(alias, "role_snowflake"):
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel.mention} "
                f" and role {alias.role_mention} deleted successfully."
            )
        else:
            msg = (
                f"Alias `{alias.alias_name}` of type "
                f"`{alias.category}` for channel {channel.mention} "
                f"deleted successfully."
            )
        await Alias.delete(**where_kwargs)
        return msg

    @classmethod
    async def create_alias(cls, alias_name, category, channel_dict, role_dict):
        where_kwargs = {}
        where_kwargs.update(channel_dict.get("columns", None))
        msg = (
            f"Alias `{alias_name}` of type `{category}` "
            f"created successfully for channel {channel_dict.get("mention", None)}."
        )
        alias = await Alias.select(category=category, **where_kwargs, singular=True)
        if alias:
            return (
                f"Alias of type `{category}` "
                f"already exists for this channel {channel_dict.get("mention", None)}."
            )
        if role_dict:
            where_kwargs.update(role_dict.get("columns", None))
            msg = (
                f"Alias `{alias_name}` of type `{category}` "
                f"created successfully for channel {channel_dict.get("mention", None)} with role {role_dict.get("mention", None)}."
            )
        alias = Alias(alias_name=alias_name, category=str(category), **where_kwargs)
        await alias.create()
        return msg
