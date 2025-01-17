```python
# pdf_manager.py
'''
pdf_manager.py

Provides PDF management functionality for the Discord bot.
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

class PDFManager:
    def __init__(self, db_pool):
        """
        :param db_pool: An asyncpg connection pool or similar.
        """
        self.pool = db_pool

    async def upload_pdf(self, reference_id: int, file_url: str) -> int:
        """
        Uploads a PDF associated with a reference.

        :return: The ID of the newly uploaded PDF.
        """
        query = """
            INSERT INTO pdfs (reference_id, file_url)
            VALUES ($1, $2)
            RETURNING id
        """
        try:
            async with self.pool.acquire() as conn:
                pdf_id = await conn.fetchval(query, reference_id, file_url)
                return pdf_id
        except Exception as e:
            logger.error(f"Failed to upload PDF: {e}")
            raise

    async def annotate_pdf(
        self,
        pdf_id: int,
        user_id: int,
        page_number: Optional[int],
        content: str,
        highlighted_text: Optional[str] = None
    ) -> int:
        """
        Adds an annotation to a PDF.

        :return: The ID of the newly created annotation.
        """
        query = """
            INSERT INTO annotations (pdf_id, user_id, page_number, content, highlighted_text)
            VALUES ($1, $2, $3, $4, $5)
            RETURNING id
        """
        try:
            async with self.pool.acquire() as conn:
                annotation_id = await conn.fetchval(
                    query, 
                    pdf_id, 
                    user_id, 
                    page_number, 
                    content, 
                    highlighted_text
                )
                return annotation_id
        except Exception as e:
            logger.error(f"Failed to annotate PDF: {e}")
            raise

    async def view_pdf(self, pdf_id: int) -> Optional[Dict]:
        """
        Retrieves information about a specific PDF.

        :return: PDF data or None if not found.
        """
        query = """
            SELECT * FROM pdfs
            WHERE id = $1
        """
        try:
            async with self.pool.acquire() as conn:
                row = await conn.fetchrow(query, pdf_id)
                return dict(row) if row else None
        except Exception as e:
            logger.error(f"Failed to view PDF: {e}")
            raise

    async def get_annotations(self, pdf_id: int) -> List[Dict]:
        """
        Retrieves all annotations for a specific PDF.

        :return: A list of annotations.
        """
        query = """
            SELECT * FROM annotations
            WHERE pdf_id = $1
        """
        try:
            async with self.pool.acquire() as conn:
                rows = await conn.fetch(query, pdf_id)
                return [dict(row) for row in rows]
        except Exception as e:
            logger.error(f"Failed to get annotations: {e}")
            raise

    async def delete_pdf(self, pdf_id: int, user_id: int) -> bool:
        """
        Deletes a PDF from the database.

        :return: True if deletion was successful, False otherwise.
        """
        query = """
            DELETE FROM pdfs
            WHERE id = $1 AND id IN (
                SELECT pdf_id FROM references 
                JOIN pdfs ON references.id = pdfs.reference_id
                WHERE references.user_id = $2
            )
        """
        try:
            async with self.pool.acquire() as conn:
                result = await conn.execute(query, pdf_id, user_id)
                return result.endswith("DELETE 1")
        except Exception as e:
            logger.error(f"Failed to delete PDF: {e}")
            raise
```