''' create_batch_completionn.py  The purpose of this program is to load a better queue for 50% cheaper API requests during low traffic times for OpenAI.
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
import json
import os
import aiohttp
import asyncio
from datetime import datetime
from lucy.utils.config import Config
from lucy.utils.helpers import *

# Load API Key
config = Config().get_config()
OPENAI_API_KEY = config["api_keys"]["OpenAI"]["api_key"]

class BatchProcessor:
    def __init__(self, bot):
        self.bot = bot  # Pass bot instance for user lookup
        self.results = {}  # Store batch results

        # Load existing results if available
        if os.path.exists(PATH_OPENAI_RESULTS):
            with open(PATH_OPENAI_RESULTS, "r") as f:
                self.results = json.load(f)

    async def upload_file(self):
        """Uploads batch request JSONL file to OpenAI using aiohttp.FormData()."""
        url = "https://api.openai.com/v1/files"
        headers = {
            "Authorization": f"Bearer {OPENAI_API_KEY}",
        }
    
        # Check if file exists
        if not os.path.exists(PATH_OPENAI_REQUESTS) or os.stat(PATH_OPENAI_REQUESTS).st_size == 0:
            return "No batch requests to process."
    
        # Prepare form data for file upload
        form_data = aiohttp.FormData()
        form_data.add_field(
            "file",
            open(PATH_OPENAI_REQUESTS, "rb"),
            filename=os.path.basename(PATH_OPENAI_REQUESTS),
            content_type="application/jsonl"
        )
        form_data.add_field("purpose", "batch")
    
        async with aiohttp.ClientSession() as session:
            async with session.post(url, headers=headers, data=form_data) as response:
                resp_json = await response.json()
                if response.status == 200:
                    return resp_json["id"]  # Return the uploaded file ID
                else:
                    return f"File upload failed: {resp_json}"

    async def create_batch(self, input_file_id):
        """Creates a batch request using OpenAI API."""
        url = "https://api.openai.com/v1/batches"
        headers = {"Authorization": f"Bearer {OPENAI_API_KEY}", "Content-Type": "application/json"}
        data = {"input_file_id": input_file_id, "endpoint": "/v1/chat/completions", "completion_window": "24h"}

        async with aiohttp.ClientSession() as session:
            async with session.post(url, headers=headers, json=data) as response:
                return await response.json()

    async def retrieve_batch_results(self, batch_id):
        """Retrieves and processes batch results."""
        url = f"https://api.openai.com/v1/batches/{batch_id}"
        headers = {"Authorization": f"Bearer {OPENAI_API_KEY}"}

        async with aiohttp.ClientSession() as session:
            while True:
                async with session.get(url, headers=headers) as response:
                    batch_status = await response.json()

                    if batch_status["status"] == "completed":
                        output_file_id = batch_status["output_file_id"]
                        return await self.download_results(output_file_id)

                    elif batch_status["status"] in ["failed", "cancelled"]:
                        return f"Batch processing failed: {batch_status}"

                    print("Batch still processing... Checking again in 5 minutes.")
                    await asyncio.sleep(300)

    async def download_results(self, output_file_id):
        """Downloads batch results and organizes them alphabetically by user."""
        url = f"https://api.openai.com/v1/files/{output_file_id}/content"
        headers = {"Authorization": f"Bearer {OPENAI_API_KEY}"}

        async with aiohttp.ClientSession() as session:
            async with session.get(url, headers=headers) as response:
                results = await response.json()

                # Organize responses alphabetically by username
                for entry in results:
                    user_id = entry["custom_id"].split('-')[0]
                    user = self.bot.get_user(int(user_id))  # Get user by ID

                    if not user:
                        continue  # Skip if user not found

                    username = user.name
                    first_letter = username[0].upper()

                    if first_letter not in self.results:
                        self.results[first_letter] = {}

                    if username not in self.results[first_letter]:
                        self.results[first_letter][username] = []

                    self.results[first_letter][username].append(entry["response"])

                # Save to file
                with open(PATH_OPENAI_RESULTS, "w") as file:
                    json.dump(self.results, file, indent=4)

                return "Batch processing completed! Results saved."

    async def process_batches(self):
        """Runs batch processing."""
        if not os.path.exists(PATH_OPENAI_REQUESTS) or os.stat(PATH_OPENAI_REQUESTS).st_size == 0:
            return "No batch requests to process."

        file_id = await self.upload_file()
        batch_response = await self.create_batch(file_id)
        batch_id = batch_response["id"]

        print(f"Batch created: {batch_id}. Waiting for results...")
        result_message = await self.retrieve_batch_results(batch_id)

        # Clear batch request file after processing
        open(PATH_OPENAI_REQUESTS, "w").close()

        return result_message

    def get_user_responses(self, user: discord.User):
        """Retrieves and deletes user's responses."""
        first_letter = user.name[0].upper()

        if first_letter in self.results and user.name in self.results[first_letter]:
            responses = self.results[first_letter][user.name]

            # Remove responses after retrieval
            del self.results[first_letter][user.name]
            if not self.results[first_letter]:  # Remove empty letters
                del self.results[first_letter]

            # Save updated results
            with open(PATH_OPENAI_RESULTS, "w") as file:
                json.dump(self.results, file, indent=4)

            return responses

        return None
