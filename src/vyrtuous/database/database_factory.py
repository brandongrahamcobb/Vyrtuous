"""stage.py The purpose of this program is to the parent class to all database operation classes.

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

from vyrtuous.bot.discord_bot import DiscordBot


class DatabaseFactory(object):

    def __init__(self):
        self.bot = DiscordBot.get_instance()

    async def create(self):
        table_name = getattr(self, "TABLE_NAME")
        fields = getattr(self, "REQUIRED_INSTANTIATION_ARGS") + getattr(
            self, "OPTIONAL_ARGS"
        )
        insert_fields = [f for f in fields if getattr(self, f, None) is not None]
        if not insert_fields:
            raise ValueError("No fields available to insert")
        placeholders = ", ".join(f"${i+1}" for i in range(len(insert_fields)))
        values = [getattr(self, f) for f in insert_fields]
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
        fields = getattr(cls, "REQUIRED_INSTANTIATION_ARGS") + getattr(
            cls, "OPTIONAL_ARGS"
        )
        table_name = getattr(cls, "TABLE_NAME")
        for key in kwargs:
            if key not in fields:
                raise ValueError(f"Invalid argument '{key}' for {cls.__name__}")
        conditions = []
        values = []
        if kwargs:
            for index, field in enumerate(sorted(kwargs)):
                conditions.append(f"{field}=${index+1}")
                values.append(kwargs[field])
        where_clause = "WHERE " + " AND ".join(conditions) if conditions else ""
        async with bot.db_pool.acquire() as conn:
            await conn.execute(f"DELETE FROM {table_name} {where_clause}", *values)

    @classmethod
    async def select(cls, *, singular=False, **kwargs):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls, "TABLE_NAME")
        fields = getattr(cls, "REQUIRED_INSTANTIATION_ARGS") + getattr(
            cls, "OPTIONAL_ARGS"
        )
        virtual_filters = {"expired"}
        for key in kwargs:
            if key not in fields and key not in virtual_filters:
                raise ValueError(f"Invalid argument '{key}' for {cls.__name__}")
        conditions = []
        values = []
        for field, value in kwargs.items():
            if field == "expired":
                if value is True:
                    conditions.append("expires_in IS NOT NULL AND expires_in < NOW()")
                continue
            conditions.append(f"{field}=${len(values)+1}")
            values.append(value)
        where_clause = "WHERE " + " AND ".join(conditions) if conditions else ""
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                f"SELECT * FROM {table_name} {where_clause}", *values
            )
        if singular:
            if not rows:
                return None
            row = rows[0]
            row_data = {k: row[k] for k in fields if k in row}
            return cls(**row_data)
        children = []
        for row in rows:
            row_data = {k: row[k] for k in fields if k in row}
            children.append(cls(**row_data))
        return children

    @classmethod
    async def update(cls, *, set_kwargs: dict, where_kwargs: dict):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls, "TABLE_NAME")
        fields = getattr(cls, "REQUIRED_INSTANTIATION_ARGS") + getattr(
            cls, "OPTIONAL_ARGS"
        )
        for key in set_kwargs.keys():
            if key not in fields:
                raise ValueError(f"Invalid update field '{key}' for {cls.__name__}")
        for key in where_kwargs.keys():
            if key not in fields:
                raise ValueError(f"Invalid filter field '{key}' for {cls.__name__}")
        if not set_kwargs:
            raise ValueError("No fields provided to update")
        if not where_kwargs:
            raise ValueError("No fields provided to filter update")
        set_fields = sorted(set_kwargs.keys())
        where_fields = sorted(where_kwargs.keys())
        assignments = [
            f"{field} = ${index + 1}" for index, field in enumerate(set_fields)
        ]
        conditions = [
            f"{field} = ${index + 1 + len(set_fields)}"
            for index, field in enumerate(where_fields)
        ]
        values = [set_kwargs[field] for field in set_fields] + [
            where_kwargs[field] for field in where_fields
        ]
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                f"""
                UPDATE {table_name}
                SET {', '.join(assignments)}
                WHERE {' AND '.join(conditions)}
            """,
                *values,
            )
