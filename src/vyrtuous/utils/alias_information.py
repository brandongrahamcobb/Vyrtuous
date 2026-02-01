from datetime import datetime, timezone

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.mgmt.alias import Alias
from vyrtuous.db.mgmt.cap import Cap
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.fields.duration import DurationObject, DurationError
from vyrtuous.utils.highest_role import resolve_highest_role


class AliasInformation:

    information = {}

    @classmethod
    async def build(cls, message):
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
        AliasInformation.information["alias"] = alias
        AliasInformation.information["snowflake_kwargs"] = {
            "channel_snowflake": int(alias.channel_snowflake),
            "guild_snowflake": int(alias.guild_snowflake),
        }
        if getattr(alias, "role_snowflake"):
            AliasInformation.information["snowflake_kwargs"].update(
                {"role_snowflake": int(alias.role_snowflake)}
            )
        AliasInformation.information["executor_role"] = await resolve_highest_role(
            channel_snowflake=int(alias.channel_snowflake),
            guild_snowflake=int(alias.guild_snowflake),
            member_snowflake=int(message.author.id),
        )
        kwargs = alias.fill_map(args)
        for field, tuple in kwargs.items():
            if field == "duration":
                duration = (
                    DurationObject(tuple[1])
                    if tuple[1] is not None
                    else DurationObject("8h")
                )
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
                        duration_str = DurationObject.from_seconds(
                            AliasInformation.information["cap_duration"]
                        )
                        AliasInformation.information["duration"] = duration_str
                        raise DurationError(information=AliasInformation.information)
                AliasInformation.information["expires_in"] = (
                    datetime.now(timezone.utc) + duration.to_timedelta()
                )
            if field == "member":
                member_dict = await do.determine_from_target(target=tuple[1])
                AliasInformation.information["snowflake_kwargs"].update(
                    {"member_snowflake": member_dict.get("id", None)}
                )
            if field == "reason":
                reason = tuple[1]
                AliasInformation.information["reason"] = reason
        return AliasInformation.information
