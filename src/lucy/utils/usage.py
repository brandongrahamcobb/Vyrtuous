import aiohttp
import asyncio
from typing import List, Optional, Any, Dict, Union
from dataclasses import dataclass, field

# Data Classes for Responses
@dataclass
class ModerationResult:
    object: str
    input_tokens: int
    num_model_requests: int
    project_id: Optional[str]
    user_id: Optional[str]
    api_key_id: Optional[str]
    model: Optional[str]

@dataclass
class ModerationBucket:
    object: str
    start_time: int
    end_time: int
    results: List[ModerationResult]

@dataclass
class CompletionResult:
    object: str
    input_tokens: int
    output_tokens: int
    input_cached_tokens: int
    input_audio_tokens: int
    output_audio_tokens: int
    num_model_requests: int
    project_id: Optional[str]
    user_id: Optional[str]
    api_key_id: Optional[str]
    model: Optional[str]
    batch: Optional[bool]

@dataclass
class CompletionBucket:
    object: str
    start_time: int
    end_time: int
    results: List[CompletionResult]

@dataclass
class UsagePage:
    object: str
    data: List[Union[ModerationBucket, CompletionBucket]]
    has_more: bool
    next_page: Optional[str]

# OpenAIUsageClient Class
class OpenAIUsageClient:
    BASE_URL = "https://api.openai.com/v1/organization/usage"

    def __init__(self, api_key: str, organization_id: Optional[str] = None):
        """
        Initialize the OpenAIUsageClient.
        :param api_key: Your OpenAI API key.
        :param organization_id: (Optional) Your OpenAI Organization ID.
        """
        self.api_key = api_key
        self.organization_id = organization_id
        self.headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }

    async def _get_usage(self, endpoint: str, params: Dict[str, Any]) -> UsagePage:
        """
        Internal method to perform GET requests to the Usage API.
        :param endpoint: The specific usage endpoint ('completions' or 'moderations').
        :param params: Query parameters as a dictionary.
        :return: Parsed UsagePage object.
        """
        url = f"{self.BASE_URL}/{endpoint}"
        if self.organization_id:
            params["organization_id"] = self.organization_id

        async with aiohttp.ClientSession() as session:
            async with session.get(url, headers=self.headers, params=params) as response:
                if response.status != 200:
                    text = await response.text()
                    raise Exception(f"API request failed with status {response.status}: {text}")
                data = await response.json()

        # Parse the response into data classes
        usage_page = UsagePage(
            object=data.get("object"),
            data=[],
            has_more=data.get("has_more"),
            next_page=data.get("next_page")
        )

        for bucket in data.get("data", []):
            if endpoint == "moderations":
                results = [ModerationResult(**result) for result in bucket.get("results", [])]
                usage_bucket = ModerationBucket(
                    object=bucket.get("object"),
                    start_time=bucket.get("start_time"),
                    end_time=bucket.get("end_time"),
                    results=results
                )
            elif endpoint == "completions":
                results = [CompletionResult(**result) for result in bucket.get("results", [])]
                usage_bucket = CompletionBucket(
                    object=bucket.get("object"),
                    start_time=bucket.get("start_time"),
                    end_time=bucket.get("end_time"),
                    results=results
                )
            else:
                continue  # Unknown endpoint
            usage_page.data.append(usage_bucket)

        return usage_page

    async def get_moderations_usage(self, **params) -> UsagePage:
        return await self._get_usage("moderations", params)

    async def get_completions_usage(self, **params) -> UsagePage:
        return await self._get_usage("completions", params)

# Example Usage
async def main():
    # Replace with your actual OpenAI API key
    api_key = "YOUR_OPENAI_API_KEY"

    # Initialize the client
    client = OpenAIUsageClient(api_key=api_key)

    # Define your query parameters
    start_time = 1730419200  # Example Unix timestamp
    limit = 1

    # Fetch Moderations usage

# To run the example, uncomment the following lines:
# if __name__ == "__main__":
#     asyncio.run(main())
