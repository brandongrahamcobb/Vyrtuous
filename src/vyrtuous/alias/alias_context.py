from datetime import datetime, timezone
from typing import Dict, Tuple

from vyrtuous.alias.alias import Alias
from vyrtuous.alias.alias_service import AliasService
# from vyrtuous.utils.permission_service import PermissionService


class AliasContext:
    MODEL = Alias

    def __init__(
        self,
        message,
        *,
        bot=None,
        cap_service=None,
        database_factory=None,
        duration=None,
    ):
        super().__init__(message=message)
        self.alias = None
        self.alias_name = None
        self.args = []
        self.expires_in = None
        self.kwargs: Dict[str, Tuple[int, str]] = {}
        self.source_kwargs: Dict[str, int] = {}
        self.message = message
        self.reason = "No reason provided"
        self.target_channel_snowflake = None
        self.target_member_snowflake = None
        self.target_role_snowflake = None
        self.record = None
        self.__bot = bot
        self.__cap_service = cap_service
        self.__database_factory = database_factory
        self.__database_factory.model = self.MODEL
        self.__duration = duration

    async def setup(self):
        self.build_source_kwargs()
        self.message_to_args()
        self.alias_name_from_args()
        await self.populate_alias()
        self.fill_map()
        await self.convert_args_to_values()

    def message_to_args(self) -> None:
        self.args = (
            self.message.content[len(self.__bot.config["discord_command_prefix"]) :]
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
        alias_entry = await self__database_factory.select(
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

    async def convert_args_to_values(self):
        for field, tuple in self.kwargs.items():
            value = tuple[1]
            if field == "duration":
                if not value:
                    duration = self.__duration("8h")
                else:
                    duration = self.__duration(value)
                # if await self.__cap_service.assert_duration_exceeds_cap(
                #     duration=duration,
                #     source_kwargs=self.source_kwargs,
                #     category=self.alias.category,
                # ):
                #     await PermissionService.check(
                #         channel_snowflake=self.target_channel_snowflake,
                #         guild_snowflake=self.source_guild_snowflake,
                #         member_snowflake=self.target_member_snowflake,
                #         lowest_role="Coordinator",
                #     )
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
