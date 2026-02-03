from datetime import datetime, timezone

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.base.alias.alias import Alias
from vyrtuous.db.base.role.role_alias import RoleAlias
from vyrtuous.db.infractions.ban.ban_alias import BanAlias
from vyrtuous.db.infractions.flag.flag_alias import FlagAlias
from vyrtuous.db.infractions.tmute.text_mute_alias import TextMuteAlias
from vyrtuous.db.infractions.vmute.voice_mute_alias import VoiceMuteAlias
from vyrtuous.db.mgmt.cap.cap import Cap
from vyrtuous.db.roles.vegan.vegan_alias import VeganAlias
from vyrtuous.fields.duration import DurationError, DurationObject
from vyrtuous.service.discord_object_service import DiscordObject
from vyrtuous.utils.check import has_equal_or_lower_role
from vyrtuous.utils.highest_role import resolve_highest_role


class AliasInformation:

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
        information = {}
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
        information["alias"] = alias_class
        information["snowflake_kwargs"] = {
            "channel_snowflake": int(alias.channel_snowflake),
            "guild_snowflake": int(alias.guild_snowflake),
            "member_snowflake": int(message.author.id),
        }
        information["executor_role"] = await resolve_highest_role(
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
                information["duration"] = duration
                cap = await Cap.select(
                    category=information["alias"].category,
                    channel_snowflake=int(alias.channel_snowflake),
                    guild_snowflake=int(alias.guild_snowflake),
                    singular=True,
                )
                if not hasattr(cap, "duration"):
                    information["cap_duration"] = DurationObject("8h").to_seconds()
                else:
                    information["cap_duration"] = cap.duration_seconds
                if (
                    duration.to_timedelta().total_seconds()
                    > information["cap_duration"]
                    or duration.number == 0
                ):
                    if information["executor_role"] == "Moderator":
                        raise DurationError(information=information)
                information["expires_in"] = (
                    None
                    if duration.number == 0
                    else datetime.now(timezone.utc) + duration.to_timedelta()
                )
            if field == "member":
                member_dict = await do.determine_from_target(target=tuple[1])
                await has_equal_or_lower_role(
                    snowflake_kwargs=information["snowflake_kwargs"],
                    member_snowflake=member_dict.get("id", None),
                )
                information["snowflake_kwargs"].update(
                    {"member_snowflake": member_dict.get("id", None)}
                )
            if field == "reason":
                reason = tuple[1]
                if not reason:
                    information["reason"] = "No reason provided."
                else:
                    information["reason"] = reason
        if getattr(alias, "role_snowflake"):
            information["snowflake_kwargs"].update(
                {"role_snowflake": int(alias.role_snowflake)}
            )
        return information
