import os

database = os.getenv("POSTGRES_DB")
host = os.getenv("POSTGRES_HOST")
password = os.getenv("POSTGRES_PASSWORD")
user = os.getenv("POSTGRES_USER")

dsn = f"postgres://{user}:{password}@{host}:{5432}/{database}"