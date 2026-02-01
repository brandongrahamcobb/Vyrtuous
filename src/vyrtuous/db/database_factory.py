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

from typing import TypeVar, Type, overload, Literal

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.utils.logger import logger

T = TypeVar("T", bound="DatabaseFactory")


class DatabaseFactory(object):

    async def create(self):
        bot = DiscordBot.get_instance()
        table_name = getattr(self, "TABLE_NAME")
        fields = getattr(self, "REQUIRED_ARGS") + getattr(
            self, "OPTIONAL_ARGS"
        )
        insert_fields = [
            f for f in fields if hasattr(self, f) and getattr(self, f) is not None
        ]
        if not insert_fields:
            raise ValueError("No fields available to insert")
        placeholders = ", ".join(f"${i+1}" for i in range(len(insert_fields)))
        values = [getattr(self, f) for f in insert_fields]
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                f"""
                INSERT INTO {table_name} ({', '.join(insert_fields)})
                VALUES ({placeholders})
                ON CONFLICT DO NOTHING
            """,
                *values,
            )
        logger.info(f"Created entry in {table_name}.")

    @classmethod
    async def delete(cls, **kwargs):
        bot = DiscordBot.get_instance()
        fields = getattr(cls, "REQUIRED_ARGS") + getattr(
            cls, "OPTIONAL_ARGS"
        )
        table_name = getattr(cls, "TABLE_NAME")
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in fields}
        conditions = []
        values = []
        if filtered_kwargs:
            for index, field in enumerate(sorted(filtered_kwargs)):
                conditions.append(f"{field}=${index+1}")
                values.append(filtered_kwargs[field])
        where_clause = "WHERE " + " AND ".join(conditions) if conditions else ""
        async with bot.db_pool.acquire() as conn:
            await conn.execute(f"DELETE FROM {table_name} {where_clause}", *values)
        logger.info(f"Deleted entry from {table_name}.")

    @classmethod
    @overload
    async def select(
        cls: Type[T], *, singular: Literal[True], inside=False, **kwargs
    ) -> T:
        ...

    @classmethod
    @overload
    async def select(
        cls: Type[T], *, singular: Literal[False], inside=False, **kwargs
    ) -> list[T]:
        ...

    @classmethod
    async def select(
        cls: Type[T], *, singular=False, inside=False, **kwargs
    ) -> T | list[T]:
        bot = DiscordBot.get_instance()
        table_name = getattr(cls, "TABLE_NAME")
        fields = getattr(cls, "REQUIRED_ARGS") + getattr(
            cls, "OPTIONAL_ARGS"
        )
        virtual_filters = {"expired"}
        real_kwargs = {k: v for k, v in kwargs.items() if k in fields}
        virtual_kwargs = {k: v for k, v in kwargs.items() if k in virtual_filters}
        conditions = []
        values = []
        if virtual_kwargs.get("expired") is True:
            conditions.append("expires_in IS NOT NULL AND expires_in < NOW()")
            real_kwargs.pop("expired", None)
        for field, value in real_kwargs.items():
            if inside:
                conditions.append(f"${len(values)+1} = ANY({field})")
            else:
                conditions.append(f"{field}=${len(values)+1}")
            values.append(value)
        where_clause = "WHERE " + " AND ".join(conditions) if conditions else ""
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                f"SELECT * FROM {table_name} {where_clause}", *values
            )
        if singular:
            if not rows:
                return []
            row = rows[0]
            row_data = {k: row[k] for k in fields if k in row}
            return cls(**row_data)
        children = []
        for row in rows:
            row_data = {k: row[k] for k in fields if k in row}
            children.append(cls(**row_data))
        logger.info(f"Selected entry from {table_name}.")
        return children

    @classmethod
    async def update(cls, *, set_kwargs: dict, where_kwargs: dict):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls, "TABLE_NAME")
        fields = getattr(cls, "REQUIRED_ARGS") + getattr(
            cls, "OPTIONAL_ARGS"
        )
        set_filtered_kwargs = {k: v for k, v in set_kwargs.items() if k in fields}
        where_filtered_kwargs = {k: v for k, v in where_kwargs.items() if k in fields}
        set_fields = sorted(set_filtered_kwargs.keys())
        where_fields = sorted(where_filtered_kwargs.keys())
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
        logger.info(f"Updated entry from {table_name}.")

    @classmethod
    async def primary_keys(cls):
        bot = DiscordBot.get_instance()
        table_name = getattr(cls, "TABLE_NAME")
        statement = """
            SELECT kcu.column_name
              FROM information_schema.table_constraints tc
              JOIN information_schema.key_column_usage kcu
                ON tc.constraint_name = kcu.constraint_name
               AND tc.table_schema = kcu.table_schema
              WHERE tc.constraint_type = 'PRIMARY KEY'
                AND tc.table_schema = 'public'
                AND tc.table_name = $1
              ORDER BY kcu.ordinal_position;
        """
        kwargs = []
        async with bot.db_pool.acquire() as conn:
            rows = await conn.fetch(statement, table_name)
        for row in rows:
            kwargs.append(row["column_name"])
        return kwargs