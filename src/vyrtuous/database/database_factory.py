from vyrtuous.bot.discord_bot import DiscordBot


class DatabaseFactory(object):

    def __init__(self):
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
        table_name = getattr(self, "TABLE_NAME")
        fields = getattr(self, "REQUIRED_INSTANTIATION_ARGS") + getattr(self, "OPTIONAL_ARGS")
        for key in kwargs:
            if key not in fields:
                raise ValueError(f"Invalid argument '{key}' for {self.__class__.__name__}")
        insert_fields = list(kwargs)
        placeholders = ", ".join(f"${i+1}" for i in range(len(insert_fields)))
        values = [kwargs[f] for f in insert_fields]
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute(
                f"""
                INSERT INTO {table_name} ({', '.join(insert_fields)})
                VALUES ({placeholders})
                ON CONFLICT DO NOTHING
                """,
                *values,
            )
    
    @classmethod
    async def delete(cls, **kwargs):
        bot = DiscordBot.get_instance()
        fields = getattr(cls, 'REQUIRED_INSTANTIATION_ARGS') + getattr(cls, 'OPTIONAL_ARGS')
        table_name = getattr(cls, 'TABLE_NAME')
        for key in kwargs:
            if key not in fields:
                raise ValueError(f"Invalid argument '{key}' for {cls.__name__}")
        conditions = []
        values = []
        if kwargs:
            for index, field in enumerate(sorted(kwargs)):
                conditions.append(f'{field}=${index+1}')
                values.append(kwargs[field])
        where_clause = 'WHERE ' + ' AND '.join(conditions) if conditions else ''
        async with bot.db_pool.acquire() as conn:
            await conn.execute(f'DELETE FROM {table_name} {where_clause}', *values)
            print(f'DELETE FROM {table_name} {where_clause}')
    
    @classmethod
    async def select(cls, **kwargs):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls, "TABLE_NAME")
        fields = getattr(cls, "REQUIRED_INSTANTIATION_ARGS") + getattr(cls, "OPTIONAL_ARGS")
        for key in kwargs:
            if key not in fields:
                raise ValueError(f"Invalid argument '{key}' for {cls.__name__}")
        conditions = []
        values = []
        for field in kwargs:
            if field == "expired" and kwargs.get("expired") is True:
                conditions.append("expires_in < NOW()")
            else:
                conditions.append(f"{field}=${len(values)+1}")
                values.append(kwargs[field])
        where_clause = "WHERE " + " AND ".join(conditions) if conditions else ''
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(f"SELECT * FROM {table_name} {where_clause}", *values)
            print(f'SELECT * FROM {table_name} {where_clause}')
        children = []
        for row in rows:
            row_data = {k: row[k] for k in row.keys() if k in fields}
            children.append(cls(**row_data))
        return children
    
    @classmethod
    async def update(cls, *, set_kwargs: dict, where_kwargs: dict):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls, 'TABLE_NAME')
        fields = getattr(cls, 'REQUIRED_INSTANTIATION_ARGS') + getattr(cls, 'OPTIONAL_ARGS')
        for key in set_kwargs.keys():
            if key not in fields:
                raise ValueError(f"Invalid update field '{key}' for {cls.__name__}")
        for key in where_kwargs.keys():
            if key not in fields:
                raise ValueError(f"Invalid filter field '{key}' for {cls.__name__}")
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
