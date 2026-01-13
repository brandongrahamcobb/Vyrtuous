from vyrtuous.bot.discord_bot import DiscordBot


class DatabaseFactory(object):

    ACT = ""
    PLURAL = ""
    SINGULAR = ""
    UNDO = ""
    REQUIRED_INSTANTIATION_ARGS = []
    OPTIONAL_ARGS = []
    TABLE_NAME = ""

    def __init(self):
        self.bot = DiscordBot.get_instance()

        # if actions:
        #     for action in actions:
        #         await History.save_entry(
        #             ctx_interaction_or_message=ctx_interaction_or_message,
        #             action_type=cls.UNDO,
        #             channel_snowflake=channel_snowflake,
        #             duration=None,
        #             highest_role=highest_role,
        #             is_modification=is_modification,
        #             member_snowflake=action.member_snowflake,
        #             reason="Clear command"
        #         )

    async def create(self, **kwargs):
        table_name = getattr(self.__class__, "TABLE_NAME")
        fields = getattr(self.__class__, "REQUIRED_INSTANTIATION_ARGS") + getattr(
            self.__class__, "OPTIONAL_ARGS"
        )
        for key in kwargs.keys():
            if key not in fields:
                raise ValueError(
                    f"Invalid argument '{key}' " f"for {self.__class__.__name__}"
                )
        insert_fields = [
            f
            for f in kwargs
            if hasattr(self.__class__, f) and getattr(self.__class__, f) is not None
        ]
        placeholders = ", ".join(f"${i+1}" for i in range(len(insert_fields)))
        values = [getattr(self.__class__, f) for f in insert_fields]
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                t"""
                    INSERT INTO {table_name} ({', '.join(insert_fields)})
                    VALUES ({placeholders})
                    ON CONFLICT DO NOTHING
                """,
                *values,
            )

    @classmethod
    async def delete(cls, **kwargs):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls.__class__, "TABLE_NAME")
        fields = getattr(cls.__class__, "REQUIRED_INSTANTIATION_ARGS") + getattr(
            cls.__class__, "OPTIONAL_ARGS"
        )
        filter_fields = [
            f
            for f in fields
            if hasattr(cls.__class__, f) and getattr(cls.__class__, f) is not None
        ]
        if not filter_fields:
            raise ValueError("No fields provided to specify deletion.")
        conditions = []
        values = []
        for i, field in enumerate(filter_fields):
            value = kwargs.get(field, getattr(cls.__class__, field))
            conditions.append(f"{field} = ${i+1}")
            values.append(value)
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                t"""
                    DELETE FROM {table_name}
                    WHERE {' AND '.join(conditions)}
                """,
                *values,
            )

    @classmethod
    async def select(cls, **kwargs):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls.__class__, "TABLE_NAME")
        fields = getattr(cls.__class__, "REQUIRED_INSTANTIATION_ARGS") + getattr(
            cls.__class__, "OPTIONAL_ARGS"
        )
        for key in kwargs.keys():
            if key not in fields:
                raise ValueError(
                    f"Invalid argument '{key}' " f"for {cls.__class__.__name__}"
                )
        filter_fields = [
            f
            for f in kwargs
            if hasattr(cls.__class__, f) and getattr(cls.__class__, f) is not None
        ]
        values = []
        if filter_fields:
            conditions = []
            for index, field in enumerate(filter_fields):
                if field == "expired" and kwargs.get("expired") is True:
                    conditions.append("expires_in < NOW()")
                else:
                    value = kwargs.get(field, getattr(cls.__class__, field))
                    conditions.append(f"{field}=${index+1}")
                    values.append(value)
            where_clause = "WHERE " + " AND ".join(conditions)
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                f"""
                SELECT * FROM {table_name}
                {where_clause}
            """,
                *values,
            )
        children = []
        for row in rows:
            row_data = {k: row[k] for k in row.keys() if k in fields}
            obj = cls.__class__(**row_data)
            children.append(obj)
        return children
    
    @classmethod
    async def update(cls, *, set_kwargs: dict, where_kwargs: dict):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls.__class__, 'TABLE_NAME')
        fields = getattr(cls.__class__, 'REQUIRED_INSTANTIATION_ARGS') + getattr(cls.__class__, 'OPTIONAL_ARGS')
        for key in set_kwargs.keys():
            if key not in fields:
                raise ValueError(f"Invalid update field '{key}' for {cls.__class__.__name__}")
        for key in where_kwargs.keys():
            if key not in fields:
                raise ValueError(f"Invalid filter field '{key}' for {cls.__class__.__name__}")
        if not set_kwargs:
            raise ValueError('No fields provided to update')
        if not where_kwargs:
            raise ValueError('No fields provided to filter update')
        set_fields = sorted(set_kwargs.keys())
        where_fields = sorted(where_kwargs.keys())
        assignments = [f"{field} = ${index + 1}" for index, field in enumerate(set_fields)]
        conditions = [f"{field} = ${index + 1 + len(set_fields)}" for index, field in enumerate(where_fields)]
        values = [set_kwargs[field] for field in set_fields] + [where_kwargs[field] for field in where_fields]
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                f"""
                UPDATE {table_name}
                SET {', '.join(assignments)}
                WHERE {' AND '.join(conditions)}
                """,
                *values,
            )
