from datetime import datetime, timezone

from vyrtuous.commands.permissions.permission_service import PermissionService
from vyrtuous.db.mgmt.cap.cap_service import CapService
from vyrtuous.commands.fields.duration import DurationObject
from vyrtuous.commands.discord_object_service import DiscordObject
from typing import Dict, Tuple
from vyrtuous.db.alias.alias_service import AliasService
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.alias.alias import Alias


class AliasContext:
    def __init__(self, message):
        self.alias = None
        self.source_channel_snowflake = message.channel.id
        self.source_guild_snowflake = message.guild.id
        self.source_member_snowflake = message.author.id
        self.alias_name = None
        self.args = []
        self.do = DiscordObject(message=message)
        self.expires_in = None
        self.kwargs: Dict[str, Tuple[int, str]] = {}
        self.source_kwargs: Dict[str, int] = {}
        self.message = message
        self.reason = "No reason provided"
        self.target_channel_snowflake = None
        self.target_member_snowflake = None
        self.target_role_snowflake = None
        self.record = None

    async def setup(self):
        self.build_source_kwargs()
        self.message_to_args()
        self.alias_name_from_args()
        await self.populate_alias()
        self.fill_map()
        await self.convert_args_to_values()

    def message_to_args(self) -> None:
        bot = DiscordBot.get_instance()
        self.args = (
            self.message.content[len(bot.config["discord_command_prefix"]) :]
            .strip()
            .split()
        )

    def fill_map(self) -> None:
        map = self.alias.ARGS_MAP
        sorted_args = sorted(map.items(), key=lambda x: x[1])
        for i, (key, pos) in enumerate(sorted_args):
            if i == len(sorted_args) - 1:
                value = (
                    " ".join(str(a) for a in self.args[pos - 1 :])
                    if len(self.args) >= pos
                    else ""
                )
            else:
                value = str(self.args[pos - 1]) if len(self.args) >= pos else ""
            self.kwargs[key] = (pos, value)

    def alias_name_from_args(self):
        self.alias_name = self.args[0]

    async def populate_alias(self):
        alias_entry = await Alias.select(
            alias_name=self.alias_name,
            guild_snowflake=self.source_guild_snowflake,
            singular=True,
        )
        if not alias_entry:
            return
        self.target_channel_snowflake = int(alias_entry.channel_snowflake)
        if getattr(alias_entry, "role_snowflake"):
            self.target_role_snowflake = int(alias_entry.role_snowflake)
        alias = AliasService.alias_category_to_alias(
            alias_category=alias_entry.category
        )
        self.alias = alias
        self.record = alias.record

    def build_source_kwargs(self):
        self.source_kwargs = {
            "channel_snowflake": self.source_channel_snowflake,
            "guild_snowflake": self.source_guild_snowflake,
            "member_snowflake": self.source_member_snowflake,
        }

    async def convert_args_to_values(self):
        for field, tuple in self.kwargs.items():
            value = tuple[1]
            if field == "duration":
                if not value:
                    duration = DurationObject("8h")
                else:
                    duration = DurationObject(value)
                if await CapService.assert_duration_exceeds_cap(
                    duration=duration,
                    source_kwargs=self.source_kwargs,
                    category=self.alias.category,
                ):
                    await PermissionService.check(
                        channel_snowflake=self.target_channel_snowflake,
                        guild_snowflake=self.source_guild_snowflake,
                        member_snowflake=self.target_member_snowflake,
                        lowest_role="Coordinator",
                    )
                self.expires_in = (
                    None
                    if duration.number == 0
                    else datetime.now(timezone.utc) + duration.to_timedelta()
                )
            elif field == "member":
                member_dict = await self.do.determine_from_target(target=value)
                self.target_member_snowflake = member_dict.get("id", None)
            elif field == "reason":
                self.reason = value
