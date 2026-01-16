import asyncpg

db_pool = asyncpg.Pool()

test_pgsql_column_type_resolution(db_pool):