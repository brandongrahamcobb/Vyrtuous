from datetime import datetime, timezone

from vyrtuous.aliases.ban_alias import BanAlias
from vyrtuous.aliases.flag_alias import FlagAlias
from vyrtuous.aliases.role_alias import RoleAlias
from vyrtuous.aliases.text_mute_alias import TextMuteAlias
from vyrtuous.aliases.vegan_alias import VeganAlias
from vyrtuous.aliases.voice_mute_alias import VoiceMuteAlias
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.fields.duration import DurationError, DurationObject
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.utils.check import has_equal_or_lower_role, check
from vyrtuous.utils.highest_role import resolve_highest_role


class AliasInformation:

    information = {}
    ALIAS_MAP = {
        "ban": BanAlias,
        "flag": FlagAlias,
        "role": RoleAlias,
        "tmute": TextMuteAlias,
        "vegan": VeganAlias,
        "vmute": VoiceMuteAlias,
    }

    @classmethod
    async def build(cls, message):
        AliasInformation.information.clear()
        bot = DiscordBot.get_instance()
        do = DiscordObject(message=message)
        args = (
            message.content[len(bot.config["discord_command_prefix"]) :].strip().split()
        )
        alias_name = args[0]
        alias = await Alias.select(
            alias_name=alias_name,
            guild_snowflake=message.guild.id,
            singular=True,
        )
        if not alias:
            return
        alias_class = cls.ALIAS_MAP[alias.category]
        kwargs = alias_class.service.fill_map(alias_class=alias_class, args=args)
        AliasInformation.information["alias"] = alias_class
        AliasInformation.information["snowflake_kwargs"] = {
            "channel_snowflake": int(alias.channel_snowflake),
            "guild_snowflake": int(alias.guild_snowflake),
            "member_snowflake": int(message.author.id),
        }
        AliasInformation.information["executor_role"] = await resolve_highest_role(
            channel_snowflake=int(alias.channel_snowflake),
            guild_snowflake=int(alias.guild_snowflake),
            member_snowflake=int(message.author.id),
        )
        for field, tuple in kwargs.items():
            if field == "duration":
                value = tuple[1]
                if not value:
                    duration = DurationObject("8h")
                else:
                    duration = DurationObject(value)
                AliasInformation.information["duration"] = duration
                cap = await Cap.select(
                    category=AliasInformation.information["alias"].category,
                    channel_snowflake=int(alias.channel_snowflake),
                    guild_snowflake=int(alias.guild_snowflake),
                    singular=True,
                )
                if not hasattr(cap, "duration"):
                    AliasInformation.information["cap_duration"] = DurationObject(
                        "8h"
                    ).to_seconds()
                else:
                    AliasInformation.information["cap_duration"] = cap.duration_seconds
                if (
                    duration.to_timedelta().total_seconds()
                    > AliasInformation.information["cap_duration"]
                    or duration.number == 0
                ):
                    if AliasInformation.information["executor_role"] == "Moderator":
                        raise DurationError(information=AliasInformation.information)
                AliasInformation.information["expires_in"] = (
                    None
                    if duration.number == 0
                    else datetime.now(timezone.utc) + duration.to_timedelta()
                )
            if field == "member":
                member_dict = await do.determine_from_target(target=tuple[1])
                await has_equal_or_lower_role(
                    snowflake_kwargs=AliasInformation.information["snowflake_kwargs"],
                    member_snowflake=member_dict.get("id", None),
                )
                AliasInformation.information["snowflake_kwargs"].update(
                    {"member_snowflake": member_dict.get("id", None)}
                )
            if field == "reason":
                reason = tuple[1]
                if not reason:
                    AliasInformation.information["reason"] = "No reason provided."
                else:
                    AliasInformation.information["reason"] = reason
        if getattr(alias, "role_snowflake"):
            AliasInformation.information["snowflake_kwargs"].update(
                {"role_snowflake": int(alias.role_snowflake)}
            )
        return AliasInformation.information
