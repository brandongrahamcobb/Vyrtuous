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
from typing import Optional, List, Dict
from .setup_logging import logger

class TagManager:
    def __init__(self, db_pool):
        self.pool = db_pool

    # -------------------------------------------------------------------------
    # 1) Add a tag (defaults to tag_type='default')
    #    If you want a looping-type tag, call with tag_type='loop'
    # -------------------------------------------------------------------------
    async def add_tag(
        self,
        name: str,
        location_id: int,
        owner_id: int,
        content: Optional[str] = None,
        attachment_url: Optional[str] = None,
        tag_type: str = 'default'  # can be 'default' or 'loop'
    ):
        if not name or not location_id or not owner_id:
            raise ValueError("Name, location_id, and owner_id are required fields.")

        # We only allow 'default' or 'loop'
        if tag_type not in ('default', 'loop'):
            raise ValueError("Invalid tag_type. Must be 'default' or 'loop'.")

        query = """
            INSERT INTO tags (name, location_id, content, attachment_url, owner_id, tag_type)
            VALUES ($1, $2, $3, $4, $5, $6)
        """
        try:
            async with self.pool.acquire() as conn:
                await conn.execute(
                    query,
                    name,
                    location_id,
                    content,
                    attachment_url,
                    owner_id,
                    tag_type
                )
        except Exception as e:
            logger.error(f"Failed to add tag: {e}")
            raise RuntimeError(
                "An error occurred while adding the tag. Ensure all fields are correct and try again."
            )

    # -------------------------------------------------------------------------
    # 2) Retrieve a tag by name & location
    #    Optionally filter by tag_type (i.e. 'default' or 'loop')
    # -------------------------------------------------------------------------
    async def get_tag(
        self,
        location_id: int,
        name: str,
        tag_type: Optional[str] = None
    ) -> Dict:
        if not location_id or not name:
            raise ValueError("Location ID and tag name are required fields.")

        # If a tag_type is given, include it in WHERE
        if tag_type:
            query = """
                SELECT * FROM tags
                WHERE location_id = $1
                  AND LOWER(name) = $2
                  AND tag_type = $3
            """
        else:
            query = """
                SELECT * FROM tags
                WHERE location_id = $1
                  AND LOWER(name) = $2
            """
        try:
            async with self.pool.acquire() as conn:
                if tag_type:
                    tag = await conn.fetchrow(query, location_id, name.lower(), tag_type)
                else:
                    tag = await conn.fetchrow(query, location_id, name.lower())
                
                if tag:
                    return dict(tag)
                raise RuntimeError(f'Tag "{name}" not found.')
        except Exception as e:
            logger.error(f"Failed to retrieve tag: {e}")
            raise RuntimeError(
                "An error occurred while retrieving the tag. Please check the location ID and name."
            )

    # -------------------------------------------------------------------------
    # 3) Update a tag (identified by name, location, owner)
    #    Note: Does not allow changing tag_type in this example.
    # -------------------------------------------------------------------------
    async def update_tag(
        self,
        name: str,
        location_id: int,
        owner_id: int,
        content: Optional[str] = None,
        attachment_url: Optional[str] = None
    ) -> int:
        if not name or not location_id or not owner_id:
            raise ValueError("Name, location_id, and owner_id are required fields.")

        query = """
            UPDATE tags
               SET content = $1,
                   attachment_url = $2
             WHERE name = $3
               AND location_id = $4
               AND owner_id = $5
        """
        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, content, attachment_url, name, location_id, owner_id)
                return int(result.split(" ")[-1])  # number of rows updated
        except Exception as e:
            logger.error(f"Failed to update tag: {e}")
            raise RuntimeError(
                "An error occurred while updating the tag. Please ensure the inputs are valid."
            )

    # -------------------------------------------------------------------------
    # 4) Delete a tag
    # -------------------------------------------------------------------------
    async def delete_tag(
        self,
        name: str,
        location_id: int,
        owner_id: int
    ) -> int:
        if not name or not location_id or not owner_id:
            raise ValueError("Name, location_id, and owner_id are required fields.")

        query = """
            DELETE FROM tags
             WHERE name = $1
               AND location_id = $2
               AND owner_id = $3
        """
        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, name, location_id, owner_id)
                return int(result.split(" ")[-1])  # number of rows deleted
        except Exception as e:
            logger.error(f"Failed to delete tag: {e}")
            raise RuntimeError(
                "An error occurred while deleting the tag. Please verify the inputs and try again."
            )

    # -------------------------------------------------------------------------
    # 5) List tags in a given guild/location
    #    Optionally filter by owner, and/or by tag_type = 'default' or 'loop'
    # -------------------------------------------------------------------------
    async def list_tags(
        self,
        location_id: int,
        owner_id: Optional[int] = None,
        tag_type: Optional[str] = None
    ) -> List[Dict]:
        if not location_id:
            raise ValueError("Location ID is a required field.")

        # Base query
        query = """SELECT name, content, attachment_url, tag_type FROM tags WHERE location_id = $1"""
        params = [location_id]

        # Filter by owner if specified
        if owner_id is not None:
            query += " AND owner_id = $2"
            params.append(owner_id)

        # Filter by tag_type if specified
        if tag_type is not None:
            if tag_type not in ('default', 'loop'):
                raise ValueError("Invalid tag_type. Must be 'default' or 'loop'.")
            if owner_id is not None:
                query += " AND tag_type = $3"
            else:
                query += " AND tag_type = $2"
            params.append(tag_type)

        try:
            async with self.pool.acquire() as conn:
                rows = await conn.fetch(query, *params)
            return [dict(r) for r in rows]
        except Exception as e:
            logger.error(f"Failed to list tags: {e}")
            raise RuntimeError(
                "An error occurred while retrieving the tag list. Please check the inputs."
            )

    # =========================================================================
    # 6) Loop configuration methods:
    #    - set_loop_config: upsert channel + enabled status in loop_configs
    #    - get_loop_config: retrieve loop config for a guild
    # =========================================================================
    async def set_loop_config(self, guild_id: int, channel_id: int, enabled: bool):
        """
        Upsert a row in the loop_configs table. If the guild row doesn't exist, create it;
        otherwise update it.
        """
        if not guild_id or not channel_id:
            raise ValueError("Guild ID and channel ID are required.")

        query = """
            INSERT INTO loop_configs (guild_id, channel_id, enabled)
            VALUES ($1, $2, $3)
            ON CONFLICT (guild_id)
            DO UPDATE SET channel_id = EXCLUDED.channel_id,
                          enabled = EXCLUDED.enabled
        """
        try:
            async with self.pool.acquire() as conn:
                await conn.execute(query, guild_id, channel_id, enabled)
        except Exception as e:
            logger.error(f"Failed to set loop config: {e}")
            raise RuntimeError(
                "An error occurred while configuring the loop. Please check the inputs."
            )

    async def get_loop_config(self, guild_id: int) -> Optional[Dict]:
        """
        Return a dict like { 'guild_id': <>, 'channel_id': <>, 'enabled': <> }
        if found, otherwise None.
        """
        if not guild_id:
            raise ValueError("Guild ID is required.")

        query = "SELECT * FROM loop_configs WHERE guild_id = $1"
        try:
            async with self.pool.acquire() as conn:
                row = await conn.fetchrow(query, guild_id)
                if row:
                    return dict(row)
                return None
        except Exception as e:
            logger.error(f"Failed to get loop config: {e}")
            raise RuntimeError(
                "An error occurred while retrieving the loop config."
            )
