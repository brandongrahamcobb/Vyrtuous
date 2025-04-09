import aiohttp
import asyncio
from typing import List, Optional, Any, Dict, Union
from dataclasses import dataclass, field

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

class OpenAIUsageClient:
    BASE_URL = 'https://api.openai.com/v1/organization/usage'

    def __init__(self, api_key: str, organization_id: Optional[str] = None):
        self.api_key = api_key
        self.organization_id = organization_id
        self.headers = {
            'Authorization': f'Bearer {self.api_key}',
            'Content-Type': 'application/json'
        }

    async def _get_usage(self, endpoint: str, params: Dict[str, Any]) -> UsagePage:
        url = f'{self.BASE_URL}/{endpoint}'
        if self.organization_id:
            params['organization_id'] = self.organization_id

        async with aiohttp.ClientSession() as session:
            async with session.get(url, headers=self.headers, params=params) as response:
                if response.status != 200:
                    text = await response.text()
                    raise Exception(f'API request failed with status {response.status}: {text}')
                data = await response.json()
        usage_page = UsagePage(
            object=data.get('object'),
            data=[],
            has_more=data.get('has_more'),
            next_page=data.get('next_page')
        )
        for bucket in data.get('data', []):
            if endpoint == 'moderations':
                results = [ModerationResult(**result) for result in bucket.get('results', [])]
                usage_bucket = ModerationBucket(
                    object=bucket.get('object'),
                    start_time=bucket.get('start_time'),
                    end_time=bucket.get('end_time'),
                    results=results
                )
            elif endpoint == 'completions':
                results = [CompletionResult(**result) for result in bucket.get('results', [])]
                usage_bucket = CompletionBucket(
                    object=bucket.get('object'),
                    start_time=bucket.get('start_time'),
                    end_time=bucket.get('end_time'),
                    results=results
                )
            else:
                continue
            usage_page.data.append(usage_bucket)
        return usage_page

    async def get_moderations_usage(self, **params) -> UsagePage:
        return await self._get_usage('moderations', params)

    async def get_completions_usage(self, **params) -> UsagePage:
        return await self._get_usage('completions', params)

async def main():
    api_key = self.config['api_keys']['OpenAI']['api_key']
    client = OpenAIUsageClient(api_key=api_key)
    start_time = 1730419200  # Example Unix timestamp
    limit = 1
    # Fetch Moderations usage

# To run the example, uncomment the following lines:
if __name__ == '__main__':
     asyncio.run(main())
