"""!/bin/python3

database_factory.py The purpose of this program is to the database factory utility class.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from typing import Any, Generic, Literal, TypeVar, overload

from vyrtuous.bot.discord_bot import DiscordBot

T = TypeVar("T", bound="DatabaseFactory")


class DatabaseFactory(Generic[T]):
    def __init__(self, model):
        self.model = model

    async def create(self, obj) -> None:
        bot: DiscordBot = DiscordBot.get_instance()
        table_name = getattr(obj.__class__, "__tablename__")
        fields = list(obj.__class__.__annotations__.keys())
        insert_fields = [
            f for f in fields if hasattr(obj, f) and getattr(obj, f) is not None
        ]
        if not insert_fields:
            raise ValueError("No fields available to insert")
        placeholders = ", ".join(f"${i + 1}" for i in range(len(insert_fields)))
        values = [getattr(obj, f) for f in insert_fields]
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                f"""
                INSERT INTO {table_name} ({", ".join(insert_fields)})
                VALUES ({placeholders})
                ON CONFLICT DO NOTHING
            """,
                *values,
            )
        bot.logger.debug(f"Created entry in {table_name}.")

    async def upsert(self, obj) -> None:
        bot: DiscordBot = DiscordBot.get_instance()
        table_name = getattr(obj.__class__, "__tablename__")
        fields = list(obj.__class__.__annotations__.keys())
        insert_fields = [
            f for f in fields if hasattr(obj, f) and getattr(obj, f) is not None
        ]
        if not insert_fields:
            raise ValueError("No fields available to insert")
        self.model = obj.__class__
        primary_keys = await self.primary_keys()
        if not primary_keys:
            raise ValueError("No primary key defined on table")
        placeholders = ", ".join(f"${i + 1}" for i in range(len(insert_fields)))
        update_fields = [f for f in insert_fields if f not in primary_keys]
        if not update_fields:
            raise ValueError("No updatable fields provided")
        update_clause = ", ".join(f"{f}=EXCLUDED.{f}" for f in update_fields)
        values = [getattr(obj, f) for f in insert_fields]
        async with bot.db_pool.acquire() as conn:
            await conn.execute(
                f"""
                INSERT INTO {table_name} ({", ".join(insert_fields)})
                VALUES ({placeholders})
                ON CONFLICT ({", ".join(primary_keys)})
                DO UPDATE SET {update_clause}
                """,
                *values,
            )
        bot.logger.debug(f"Upserted entry in {table_name}.")

    async def delete(self, **kwargs) -> None:
        bot: DiscordBot = DiscordBot.get_instance()
        fields = list(self.model.__annotations__.keys())
        table_name = getattr(self.model, "__tablename__")
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in fields}
        conditions: list[str] = []
        values: list[object] = []
        for field in sorted(filtered_kwargs):
            value = filtered_kwargs[field]
            if value is None:
                conditions.append(f"{field} IS NULL")
            else:
                values.append(value)
                conditions.append(f"{field}=${len(values)}")
        where_clause = "WHERE " + " AND ".join(conditions) if conditions else ""
        async with bot.db_pool.acquire() as conn:
            await conn.execute(f"DELETE FROM {table_name} {where_clause}", *values)
        bot.logger.debug(f"Deleted entry from {table_name}.")

    async def delete_by_cls(self, cls, **kwargs) -> None:
        bot: DiscordBot = DiscordBot.get_instance()
        fields = list(cls.__annotations__.keys())
        table_name = getattr(cls, "__tablename__")
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in fields}
        if not filtered_kwargs:
            raise ValueError(
                f"delete_by_cls called with no matching fields for {table_name}; refusing unbounded DELETE"
            )
        conditions = []
        values = []
        if filtered_kwargs:
            for index, field in enumerate(sorted(filtered_kwargs)):
                conditions.append(f"{field}=${index + 1}")
                values.append(filtered_kwargs[field])
        where_clause = "WHERE " + " AND ".join(conditions) if conditions else ""
        async with bot.db_pool.acquire() as conn:
            await conn.execute(f"DELETE FROM {table_name} {where_clause}", *values)
        bot.logger.debug(f"Deleted entry from {table_name}.")

    @overload
    async def select(
        self, *, singular: Literal[True], inside_fields=[], **kwargs
    ) -> T: ...

    @overload
    async def select(
        self, *, singular: Literal[False], inside_fields=[], **kwargs
    ) -> list[T]: ...

    async def select(
        self, *, singular=False, inside_fields=[], **kwargs
    ) -> T | list[T]:
        bot: DiscordBot = DiscordBot.get_instance()
        table_name = getattr(self.model, "__tablename__")
        fields = list(self.model.__annotations__.keys())
        virtual_filters = {"expired"}
        real_kwargs = {k: v for k, v in kwargs.items() if k in fields}
        virtual_kwargs = {k: v for k, v in kwargs.items() if k in virtual_filters}
        conditions: list[str] = []
        values: list[Any] = []
        if virtual_kwargs.get("expired") is True:
            conditions.append("expires_in IS NOT NULL AND expires_in < NOW()")
            real_kwargs.pop("expired", None)
        for field, value in real_kwargs.items():
            if value is None:
                conditions.append(f"{field} IS NULL")
            elif field in inside_fields:
                values.append(value)
                conditions.append(f"${len(values)} = ANY({field})")
            else:
                values.append(value)
                conditions.append(f"{field}=${len(values)}")
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
            return self.model(**row_data)
        children = []
        for row in rows:
            row_data = {k: row[k] for k in fields if k in row}
            children.append(self.model(**row_data))
        bot.logger.debug(f"Selected entry from {table_name}.")
        return children

    async def update(self, *, set_kwargs: dict, where_kwargs: dict) -> None:
        bot: DiscordBot = DiscordBot.get_instance()
        table_name = getattr(self.model, "__tablename__")
        fields = list(self.model.__annotations__.keys())
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
                SET {", ".join(assignments)}
                WHERE {" AND ".join(conditions)}
            """,
                *values,
            )
        bot.logger.debug(f"Updated entry from {table_name}.")

    async def primary_keys(self) -> list[str]:
        bot: DiscordBot = DiscordBot.get_instance()
        table_name = getattr(self.model, "__tablename__")
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
