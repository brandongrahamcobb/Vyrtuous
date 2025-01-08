''' tag.py  The purpose of this program is to provide tagging and database functionality from cd ../.
    Copyright (C) 2024  github.com/brandongrahamcobb

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
'''
# tag.py
from typing import Optional, List, Dict
from .setup_logging import logger  # Or use `import logging` and `logging.getLogger(__name__)`

class TagManager:
    def __init__(self, db_pool):
        """
        :param db_pool: An asyncpg connection pool or similar.
        """
        self.pool = db_pool

    async def add_tag(
        self, 
        name: str, 
        location_id: int, 
        owner_id: int, 
        content: Optional[str] = None, 
        attachment_url: Optional[str] = None, 
        tag_type: str = 'default'
    ):
        """
        Inserts a new tag if it doesn't already exist in this location_id.
        """
        query_check = """
            SELECT 1 
            FROM tags 
            WHERE name = $1 AND location_id = $2 AND owner_id = $3
        """
        query_insert = """
            INSERT INTO tags (name, location_id, content, attachment_url, owner_id, tag_type)
            VALUES ($1, $2, $3, $4, $5, $6)
        """
        try:
            async with self.pool.acquire() as conn:
                # Check for duplicate tag
                existing_tag = await conn.fetchval(query_check, name, location_id, owner_id)
                if existing_tag:
                    raise ValueError(f"A tag with the name '{name}' already exists for you in this location.")

                # Insert the tag
                await conn.execute(query_insert, name, location_id, content, attachment_url, owner_id, tag_type)

        except Exception as e:
            logger.error(f"Failed to add tag: {e}")
            raise

    async def get_tag(self, location_id: int, name: str):
        """
        Returns a tag dict {name, content, attachment_url, tag_type} or None.
        """
        query = """
            SELECT name, content, attachment_url, tag_type
            FROM tags
            WHERE location_id = $1 
              AND LOWER(name) = LOWER($2)
        """
        try:
            async with self.pool.acquire() as conn:
                row = await conn.fetchrow(query, location_id, name)
                return dict(row) if row else None
        except Exception as e:
            logger.error(f"Failed to fetch tag: {e}")
            raise

    async def rename_tag(self, old_name: str, new_name: str, location_id: int, owner_id: int) -> int:
        """
        Rename an existing tag from old_name to new_name.
        Returns the number of rows updated (0 if none).
        Raises ValueError if the new_name already exists or an error occurs.
        """
        # 1) Check if a tag with new_name already exists in this location for the same owner
        query_check = """
            SELECT 1 
            FROM tags
            WHERE lower(name) = lower($1) 
              AND location_id = $2
              AND owner_id = $3
        """
        
        # 2) Rename query
        query_rename = """
            UPDATE tags
            SET name = $1
            WHERE lower(name) = lower($2)
              AND location_id = $3
              AND owner_id = $4
        """
        try:
            async with self.pool.acquire() as conn:
                # Check for duplicate
                existing = await conn.fetchval(query_check, new_name, location_id, owner_id)
                if existing:
                    raise ValueError(f"A tag named '{new_name}' already exists for you in this location.")
                
                # Perform the rename
                result = await conn.execute(query_rename, new_name, old_name, location_id, owner_id)
                row_count = int(result.split()[-1])  # e.g. 'UPDATE 1' → 1
                return row_count
        except Exception as e:
            logger.error(f"Failed to rename tag: {e}")
            raise

    async def update_tag(
        self,
        name: str,
        location_id: int,
        owner_id: int,
        updates: dict
    ) -> int:
        """
        Updates fields in a tag. Returns the number of rows updated (0 if not found).
        `updates` can contain keys in [content, attachment_url, tag_type].
        """
        fields_to_update = []
        params = []
        i = 1

        for field, value in updates.items():
            fields_to_update.append(f"{field} = ${i}")
            params.append(value)
            i += 1

        query = f"""
            UPDATE tags
            SET {', '.join(fields_to_update)}
            WHERE name = ${i}
              AND location_id = ${i+1}
              AND owner_id = ${i+2}
        """
        params.extend([name, location_id, owner_id])

        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, *params)
                # result is something like 'UPDATE <count>', so we can parse rowcount
                row_count = int(result.split()[-1])
                return row_count
        except Exception as e:
            logger.error(f"Failed to update tag: {e}")
            raise RuntimeError("An error occurred while updating the tag.")

    async def delete_tag(self, name: str, location_id: int, owner_id: int) -> int:
        """
        Deletes the tag from the DB. Returns the number of rows deleted (0 if none).
        """
        query = """
            DELETE FROM tags 
            WHERE name = $1 
              AND location_id = $2 
              AND owner_id = $3
        """
        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, name, location_id, owner_id)
                row_count = int(result.split()[-1])
                return row_count
        except Exception as e:
            logger.error(f"Failed to delete tag: {e}")
            raise

    async def list_tags(self, location_id: int, owner_id: Optional[int] = None, tag_type: Optional[str] = None) -> List[dict]:
        """
        Returns a list of dicts for tags that match the filters.
        If `owner_id` is provided, returns only that owner's tags.
        If `tag_type` is provided, filters by 'default' or 'loop'.
        """
        query = "SELECT name, content, attachment_url, tag_type FROM tags WHERE location_id = $1"
        params = [location_id]
        idx = 2

        if owner_id is not None:
            query += f" AND owner_id = ${idx}"
            params.append(owner_id)
            idx += 1

        if tag_type is not None:
            query += f" AND tag_type = ${idx}"
            params.append(tag_type)
            idx += 1

        try:
            async with self.pool.acquire() as conn:
                rows = await conn.fetch(query, *params)
                return [dict(row) for row in rows]
        except Exception as e:
            logger.error(f"Failed to list tags: {e}")
            raise

    async def set_loop_config(self, guild_id: int, channel_id: Optional[int], enabled: bool):
        """
        Inserts or updates a record for loop-config in the loop_configs table.
        """
        query = """
            INSERT INTO loop_configs (guild_id, channel_id, enabled)
            VALUES ($1, $2, $3)
            ON CONFLICT (guild_id)
            DO UPDATE SET 
                channel_id = EXCLUDED.channel_id,
                enabled = EXCLUDED.enabled
        """
        try:
            async with self.pool.acquire() as conn:
                await conn.execute(query, guild_id, channel_id, enabled)
        except Exception as e:
            logger.error(f"Failed to set loop config: {e}")
            raise

    async def get_loop_config(self, guild_id: int) -> Optional[dict]:
        """
        Returns a dict {guild_id, channel_id, enabled} or None if not set.
        """
        query = "SELECT guild_id, channel_id, enabled FROM loop_configs WHERE guild_id = $1"
        try:
            async with self.pool.acquire() as conn:
                row = await conn.fetchrow(query, guild_id)
                return dict(row) if row else None
        except Exception as e:
            logger.error(f"Failed to get loop config: {e}")
            raise
