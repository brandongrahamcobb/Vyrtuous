'''
reference_manager.py Provides reference management functionality for the Discord bot.
Copyright (C) 2024 github.com/brandongrahamcobb

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

from typing import List, Optional, Dict
from lucy.utils.setup_logging import logger

class ReferenceManager:
    def __init__(self, db_pool):
        """
        :param db_pool: An asyncpg connection pool or similar.
        """
        self.pool = db_pool

    async def add_reference(
        self,
        user_id: int,
        location_id: int,
        title: str,
        authors: List[str],
        publication_year: Optional[int] = None,
        doi: Optional[str] = None,
        abstract: Optional[str] = None,
        tags: Optional[List[int]] = None
    ) -> int:
        """
        Adds a new reference to the database.

        :return: The ID of the newly created reference.
        """
        query = """
            INSERT INTO reference_list (user_id, location_id, title, authors, publication_year, doi, abstract, tags)
            VALUES ($1, $2, $3, $4, $5, $6, $7, $8)
            RETURNING id
        """
        try:
            async with self.pool.acquire() as conn:
                return await conn.fetchval(
                    query, user_id, location_id, title, authors, publication_year, doi, abstract, tags
                )
        except Exception as e:
            logger.error(f"Failed to add reference (title: {title}, user_id: {user_id}): {e}")
            raise RuntimeError("Failed to add reference.") from e

    async def delete_reference(self, reference_id: int, user_id: int) -> bool:
        """
        Deletes a reference from the database.

        :return: True if deletion was successful, False otherwise.
        """
        query = "DELETE FROM reference_list WHERE id = $1 AND user_id = $2"
        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, reference_id, user_id)
                return result.endswith("DELETE 1")
        except Exception as e:
            logger.error(f"Failed to delete reference (id: {reference_id}): {e}")
            raise RuntimeError("Failed to delete reference.") from e

    async def list_references(
        self, 
        user_id: int, 
        location_id: int, 
        tag_ids: Optional[List[int]] = None
    ) -> List[Dict]:
        """
        Lists all references for a user, optionally filtered by tags.

        :return: A list of references.
        """
        base_query = "SELECT * FROM reference_list WHERE user_id = $1 AND location_id = $2"
        params = [user_id, location_id]
        if tag_ids:
            base_query += " AND tags && $3"
            params.append(tag_ids)
        try:
            async with self.pool.acquire() as conn:
                rows = await conn.fetch(base_query, *params)
                return [dict(row) for row in rows]
        except Exception as e:
            logger.error(f"Failed to list references (user_id: {user_id}): {e}")
            raise RuntimeError("Failed to list references.") from e

    async def search_references(
        self, 
        user_id: int, 
        location_id: int, 
        query_text: str
    ) -> List[Dict]:
        """
        Searches references based on a query string matching title or authors.

        :return: A list of matching references.
        """
        search_query = """
            SELECT * FROM reference_list
            WHERE user_id = $1 AND location_id = $2
              AND (
                  title ILIKE '%' || $3 || '%' 
                  OR EXISTS (
                      SELECT 1 FROM unnest(authors) AS author WHERE author ILIKE '%' || $3 || '%'
                  )
              )
        """
        try:
            async with self.pool.acquire() as conn:
                rows = await conn.fetch(search_query, user_id, location_id, query_text)
                return [dict(row) for row in rows]
        except Exception as e:
            logger.error(f"Failed to search references (query: {query_text}): {e}")
            raise RuntimeError("Reference search failed.") from e
