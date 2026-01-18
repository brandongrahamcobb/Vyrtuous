from datetime import date, datetime, timedelta
from uuid import UUID
import asyncio
import time

import asyncpg
from pglast import parse_sql
from pglast.ast import ColumnDef, Constraint, CreateStmt

from vyrtuous.database.database import Database


class NoDataFoundError(asyncpg.NoDataFoundError):

    def __init__(self, statement: str):
        self.statement = statement
        super().__init__(f"No rows found for statement: {self.statement}.")


async def main():
    db_pool = await Database().database_init()

    select_table_names_statement = """
        SELECT
            table_schema,
            table_name
        FROM information_schema.tables
        WHERE table_type = 'BASE TABLE'
          AND table_schema NOT IN ('pg_catalog', 'information_schema')
        ORDER BY table_schema, table_name;
    """

    select_table_attributes_statement = """
        SELECT
            c.column_default,
            c.column_name,
            c.data_type,
            c.is_nullable,
            fk.foreign_table,
            fk.foreign_column,
            EXISTS (
                SELECT 1
                FROM information_schema.table_constraints tc
                JOIN information_schema.key_column_usage kcu
                  ON tc.constraint_name = kcu.constraint_name
                 AND tc.table_schema = kcu.table_schema
                WHERE tc.constraint_type = 'PRIMARY KEY'
                  AND tc.table_schema = 'public'
                  AND tc.table_name = c.table_name
                  AND kcu.column_name = c.column_name
            ) AS is_primary_key,
            (
                SELECT tc.constraint_name
                FROM information_schema.table_constraints tc
                JOIN information_schema.key_column_usage kcu
                  ON tc.constraint_name = kcu.constraint_name
                 AND tc.table_schema = kcu.table_schema
                WHERE tc.constraint_type = 'PRIMARY KEY'
                  AND tc.table_schema = 'public'
                  AND tc.table_name = c.table_name
                  AND kcu.column_name = c.column_name
                LIMIT 1
            ) AS primary_key_constraint_name
        FROM information_schema.columns c
        LEFT JOIN (
            SELECT
                kcu.table_name,
                kcu.column_name,
                ccu.table_name AS foreign_table,
                ccu.column_name AS foreign_column
            FROM information_schema.table_constraints tc
            JOIN information_schema.key_column_usage kcu
                ON tc.constraint_name = kcu.constraint_name
               AND tc.table_schema = kcu.table_schema
            JOIN information_schema.constraint_column_usage ccu
                ON ccu.constraint_name = tc.constraint_name
               AND ccu.table_schema = tc.table_schema
            WHERE tc.constraint_type = 'FOREIGN KEY'
              AND tc.table_schema = 'public'
        ) fk
            ON fk.table_name = c.table_name
           AND fk.column_name = c.column_name
        WHERE c.table_schema = 'public'
          AND c.table_name = $1
        ORDER BY c.ordinal_position;
    """

    async def test_main():
        try:
            ttn_list = await test_query_table_names(
                db_pool=db_pool, statement=select_table_names_statement
            )
            tsql_dict = await test_pgsql_column_type_resolution(
                db_pool=db_pool,
                table_names=ttn_list,
                statement=select_table_attributes_statement,
            )
            for table_name, table_dict in tsql_dict.items():
                pk_columns = set()
                for column_name, column_data in table_dict.items():
                    if column_data["is_primary_key"]:
                        pk_columns.add(column_name)
                print(
                    f"""
                    TABLE: {table_name}\n
                    COLUMNS: {', '.join(tsql_dict[table_name].keys())}\n
                    PRIMARY_KEYS: ({', '.join(pk_columns)})
                """
                )
            psql_dict = test_parse_sql()
            assert psql_dict == tsql_dict
        except Exception as e:
            return print(str(e).capitalize())

    await test_main()


def test_parse_sql():
    with open("tests/unit_tests/0001-init.sql", "r") as f:
        sql_text = f.read()
        sql_dict = {}
        tree = parse_sql(sql_text)
        for node in tree:
            stmt = node.stmt
            if not isinstance(stmt, CreateStmt):
                continue
            tn = stmt.relation.relname
            sql_dict.setdefault(tn, {})
            pk_columns = set()
            for elt in stmt.tableElts:
                if isinstance(elt, Constraint) and elt.contype == "CONSTR_PRIMARY":
                    for key in elt.keys:
                        pk_columns.add(key.val)
            for elt in stmt.tableElts:
                if not isinstance(elt, ColumnDef):
                    continue
                cn = elt.colname
                sql_dict[tn].setdefault(cn, {})
                dt = " ".join(
                    (p.val if hasattr(p, "val") else str(p)) for p in elt.typeName.names
                )
                match dt:
                    case "smallint" | "integer" | "bigint":
                        sql_dict[tn][cn][dt] = int
                    case "real" | "double precision" | "numeric":
                        sql_dict[tn][cn][dt] = float
                    case "character varying" | "character" | "text":
                        sql_dict[tn][cn][dt] = str
                    case "boolean":
                        sql_dict[tn][cn][dt] = bool
                    case "date":
                        sql_dict[tn][cn][dt] = date
                    case "timestamp without time zone" | "timestamp with time zone":
                        sql_dict[tn][cn][dt] = datetime
                    case "time without time zone" | "time with time zone":
                        sql_dict[tn][cn][dt] = time
                    case "interval":
                        sql_dict[tn][cn][dt] = timedelta
                    case "uuid":
                        sql_dict[tn][cn][dt] = UUID
                    case "json" | "jsonb":
                        sql_dict[tn][cn][dt] = dict
                    case "bytea":
                        sql_dict[tn][cn][dt] = bytes
                    case _:
                        sql_dict[tn][cn][dt] = object
                sql_dict[tn][cn]["column_default"] = elt.raw_default
                sql_dict[tn][cn]["is_primary_key"] = cn in pk_columns or any(
                    c.contype == "CONSTR_PRIMARY" for c in elt.constraints or []
                )
                sql_dict[tn][cn]["is_nullable"] = not any(
                    c.contype == "CONSTR_NOTNULL" for c in elt.constraints or []
                )
        return sql_dict


async def test_pgsql_column_type_resolution(
    db_pool: asyncpg.Pool, statement: str, table_names: list[str]
):
    sql_dict = {}
    for tn in table_names:
        try:
            async with db_pool.acquire() as conn:
                rows = await conn.fetch(statement, tn)
        except Exception as e:
            raise e
        if not rows:
            raise NoDataFoundError()

        sql_dict.setdefault(tn, {})

        for row in rows:

            cn = row["column_name"]
            sql_dict[tn].setdefault(cn, {})

            cd = row["column_default"]
            if cd is None:
                sql_dict[tn][cn]["column_default"] = None
            elif cd.startswith("'") and "'::" in cd:
                sql_dict[tn][cn]["column_default"] = cd.split("'")[1]
            else:
                sql_dict[tn][cn]["column_default"] = cd

            dt = row["data_type"]
            match dt:
                case "smallint" | "integer" | "bigint":
                    sql_dict[tn][cn][dt] = int
                case "real" | "double precision" | "numeric":
                    sql_dict[tn][cn][dt] = float
                case "character varying" | "character" | "text":
                    sql_dict[tn][cn][dt] = str
                case "boolean":
                    sql_dict[tn][cn][dt] = bool
                case "date":
                    sql_dict[tn][cn][dt] = date
                case "timestamp without time zone" | "timestamp with time zone":
                    sql_dict[tn][cn][dt] = datetime
                case "time without time zone" | "time with time zone":
                    sql_dict[tn][cn][dt] = time
                case "interval":
                    sql_dict[tn][cn][dt] = timedelta
                case "uuid":
                    sql_dict[tn][cn][dt] = UUID
                case "json" | "jsonb":
                    sql_dict[tn][cn][dt] = dict
                case "bytea":
                    sql_dict[tn][cn][dt] = bytes
                case "ARRAY":
                    sql_dict[tn][cn][dt] = list
                case "USER-DEFINED":
                    sql_dict[tn][cn][dt] = str
                case _:
                    sql_dict[tn][cn][dt] = object

            is_pk = row["is_primary_key"]
            sql_dict[tn][cn]["is_primary_key"] = is_pk
            is_nullable = row["is_nullable"]
            sql_dict[tn][cn]["is_nullable"] = bool(is_nullable)
    return sql_dict


async def query_table_names(db_pool: asyncpg.Pool, statement: str) -> list[str]:
    table_names = set()
    try:
        async with db_pool.acquire() as conn:
            rows = await conn.fetch(statement)
    except Exception as e:
        raise e
    if rows:
        for row in rows:
            table_names.add(row["table_name"])
    else:
        raise NoDataFoundError()
    return sorted(list(table_names))


async def test_query_table_names(db_pool: asyncpg.Pool, statement: str):
    try:
        tn_list = await query_table_names(db_pool=db_pool, statement=statement)
        assert tn_list is not None
        return tn_list
    except Exception as e:
        raise e


if __name__ == "__main__":
    asyncio.run(main())
