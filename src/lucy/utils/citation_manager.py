```python
# citation_manager.py
'''
citation_manager.py

Provides citation management functionality for the Discord bot.
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

from typing import Optional, List, Dict
from lucy.utils.setup_logging import logger

class CitationManager:
    def __init__(self, db_pool):
        """
        :param db_pool: An asyncpg connection pool or similar.
        """
        self.pool = db_pool

    async def generate_citation(
        self,
        reference_id: int,
        user_id: int,
        citation_style: str
    ) -> int:
        """
        Generates a citation in the specified style for a reference.

        :return: The ID of the newly created citation.
        """
        # Placeholder for actual citation generation logic
        # You might integrate with a library like CiteProc-py or similar
        citation_text = f"Generated citation for reference {reference_id} in style {citation_style}"

        query = """
            INSERT INTO citations (reference_id, user_id, citation_style, citation_text)
            VALUES ($1, $2, $3, $4)
            RETURNING id
        """
        try:
            async with self.pool.acquire() as conn:
                citation_id = await conn.fetchval(
                    query, 
                    reference_id, 
                    user_id, 
                    citation_style, 
                    citation_text
                )
                return citation_id
        except Exception as e:
            logger.error(f"Failed to generate citation: {e}")
            raise

    async def export_bibliography(
        self,
        user_id: int,
        location_id: int,
        citation_style: str
    ) -> str:
        """
        Exports the user's bibliography in the specified citation style.

        :return: The formatted bibliography as a string.
        """
        query = """
            SELECT citations.citation_text
            FROM citations
            JOIN references ON citations.reference_id = references.id
            WHERE citations.user_id = $1 
              AND references.location_id = $2
              AND citations.citation_style = $3
            ORDER BY citations.created_at ASC
        """
        try:
            async with self.pool.acquire() as conn:
                rows = await conn.fetch(query, user_id, location_id, citation_style)
                bibliography = "\n".join([row['citation_text'] for row in rows])
                return bibliography
        except Exception as e:
            logger.error(f"Failed to export bibliography: {e}")
            raise

    async def delete_citation(self, citation_id: int, user_id: int) -> bool:
        """
        Deletes a citation from the database.

        :return: True if deletion was successful, False otherwise.
        """
        query = """
            DELETE FROM citations
            WHERE id = $1 AND user_id = $2
        """
        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, citation_id, user_id)
                return result.endswith("DELETE 1")
        except Exception as e:
            logger.error(f"Failed to delete citation: {e}")
            raise
```