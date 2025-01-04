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

    async def add_tag(
        self,
        name: str,
        location_id: int,
        owner_id: int,
        content: Optional[str] = None,
        attachment_url: Optional[str] = None,
    ):
        if not name or not location_id or not owner_id:
            raise ValueError("Name, location_id, and owner_id are required fields.")

        query = """
            INSERT INTO tags (name, location_id, content, attachment_url, owner_id)
            VALUES ($1, $2, $3, $4, $5)
        """
        try:
            async with self.pool.acquire() as conn:
                await conn.execute(query, name, location_id, content, attachment_url, owner_id)
        except Exception as e:
            logger.error(f"Failed to add tag: {e}")
            raise RuntimeError("An error occurred while adding the tag. Ensure all fields are correct and try again.")

    async def get_tag(self, location_id: int, name: str) -> Dict:
        if not location_id or not name:
            raise ValueError("Location ID and tag name are required fields.")
        query = """SELECT * FROM tags WHERE location_id = $1 AND LOWER(name) = $2"""
        try:
            async with self.pool.acquire() as conn:
                tag = await conn.fetchrow(query, location_id, name.lower())
                if tag:
                    return dict(tag)
                raise RuntimeError(f'Tag "{name}" not found.')
        except Exception as e:
            logger.error(f"Failed to retrieve tag: {e}")
            raise RuntimeError("An error occurred while retrieving the tag. Please check the location ID and name.")

    async def update_tag(
        self,
        name: str,
        location_id: int,
        owner_id: int,
        content: Optional[str] = None,
        attachment_url: Optional[str] = None,
    ) -> int:
        if not name or not location_id or not owner_id:
            raise ValueError("Name, location_id, and owner_id are required fields.")
        query = """
            UPDATE tags
            SET content = $1, attachment_url = $2
            WHERE name = $3 AND location_id = $4 AND owner_id = $5
        """
        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, content, attachment_url, name, location_id, owner_id)
                return int(result.split(" ")[-1])
        except Exception as e:
            logger.error(f"Failed to update tag: {e}")
            raise RuntimeError("An error occurred while updating the tag. Please ensure the inputs are valid.")

    async def delete_tag(self, name: str, location_id: int, owner_id: int) -> int:
        if not name or not location_id or not owner_id:
            raise ValueError("Name, location_id, and owner_id are required fields.")
        query = """DELETE FROM tags WHERE name = $1 AND location_id = $2 AND owner_id = $3"""
        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, name, location_id, owner_id)
                return int(result.split(" ")[-1])
        except Exception as e:
            logger.error(f"Failed to delete tag: {e}")
            raise RuntimeError("An error occurred while deleting the tag. Please verify the inputs and try again.")

    async def list_tags(self, location_id: int, owner_id: int) -> List[Dict]:
        if not location_id or not owner_id:
            raise ValueError("Location ID and owner ID are required fields.")
        query = """SELECT name, content, attachment_url FROM tags WHERE location_id = $1 AND owner_id = $2"""
        try:
            async with self.pool.acquire() as conn:
                tags = await conn.fetch(query, location_id, owner_id)
            return [dict(tag) for tag in tags]
        except Exception as e:
            logger.error(f"Failed to list tags: {e}")
            raise RuntimeError("An error occurred while retrieving the tag list. Please check the inputs.")

