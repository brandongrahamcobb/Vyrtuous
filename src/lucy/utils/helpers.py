''' helpers.py  The purpose of this program is to provide generic parameters.
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

from os.path import dirname, abspath, expanduser, join
from utils.setup_logging import logger

# Base and Home Paths
DIR_BASE = dirname(abspath(__file__))
PATH_HOME = expanduser('~')

# Config and Preferences
PATH_CONFIG_YAML = join(PATH_HOME, '.config', 'spawd', 'config.yaml')
PATH_CONFIG_JSON = join(PATH_HOME, '.config', 'spawd', 'config.json')

# Logs
PATH_LOG = join(PATH_HOME, '.log', 'spawd', 'discord.log')

# Script Paths
PATH_ADD_WATERMARK = join(DIR_BASE, 'add_watermark.py')
PATH_ADJUST_HUE_AND_SATURATION = join(DIR_BASE, 'adjust_hue_and_saturation.py')
PATH_ARPP = join(DIR_BASE, 'api_request_parallel_processor.py')
PATH_BENCHMARK = join(DIR_BASE, 'benchmark.py')
PATH_CLEAR_SCREEN = join(DIR_BASE, 'clear_screen.py')
PATH_COMBINE = join(DIR_BASE, 'combine.py')
PATH_CONFIG = join(DIR_BASE, 'config.py')
PATH_CREATE_BATCH_COMPLETION = join(DIR_BASE, 'create_batch_completion.py')
PATH_CREATE_COMPLETION_DEPRECATED = join(DIR_BASE, 'create_completion_deprecated.py')
PATH_CREATE_COMPLETION = join(DIR_BASE, 'create_completion.py')
PATH_CREATE_HTTPS_COMPLETION = join(DIR_BASE, 'create_https_completion.py')
PATH_CREATE_HTTPS_MODERATION = join(DIR_BASE, 'create_https_moderation.py')
PATH_CREATE_MODERATION = join(DIR_BASE, 'create_moderation.py')
PATH_DISCORD = join(DIR_BASE, 'discord.py')
PATH_DRAW_FINGERPRINT = join(DIR_BASE, 'draw_fingerprint.py')
PATH_DRAW_WATERMARKED_MOLECULE = join(DIR_BASE, 'draw_watermarked_molecule.py')
PATH_FINE_TUNING = join(DIR_BASE, 'fine_tuning.py')
PATH_FORMAT_ERROR_CHECK = join(DIR_BASE, 'format_error_check.py')
PATH_GET_MOLECULE_NAME = join(DIR_BASE, 'get_molecule_name.py')
PATH_GET_MOL = join(DIR_BASE, 'get_mol.py')
PATH_GET_PROXIMITY = join(DIR_BASE, 'get_proximity.py')
PATH_GET_SCRIPTURE = join(DIR_BASE, 'get_scripture.py')
PATH_GOOGLE = join(DIR_BASE, 'google.py')
PATH_GSRS = join(DIR_BASE, 'gsrs.py')
PATH_HELPERS = join(DIR_BASE, 'helpers.py')
PATH_HYBRID = join(DIR_BASE, '..', 'bot', 'cogs', 'hybrid.py')
PATH_INCREMENT_VERSION = join(DIR_BASE, 'increment_version.py')
PATH_INDICA = join(DIR_BASE, '..', 'bot', 'cogs', 'indica.py')
PATH_LINKEDIN = join(DIR_BASE, 'linkedin.py')
PATH_LOAD_CONTENTS = join(DIR_BASE, 'load_contents.py')
PATH_LOAD_YAML = join(DIR_BASE, 'load_yaml.py')
PATH_NLP_UTILS = join(DIR_BASE, 'nlp_utils.py')
PATH_PROMPT_FOR_VALUES = join(DIR_BASE, 'prompt_for_values.py')
PATH_SATIVA = join(DIR_BASE, '..', 'bot', 'cogs', 'sativa.py')
PATH_SCRIPT = join(DIR_BASE, 'script.py')
PATH_SETUP_LOGGING = join(DIR_BASE, 'setup_logging.py')
PATH_TAG = join(DIR_BASE, 'tag.py')
PATH_TRAINING = join(DIR_BASE, '..', 'training.jsonl')
PATH_UNIQUE_PAIRS = join(DIR_BASE, 'unique_pairs.py')

# Discord
DISCORD_CHARACTER_LIMITS = [2000, 4000]
DISCORD_CHARACTER_LIMIT = 2000
DISCORD_COGS = [
    'bot.cogs.hybrid',
    'bot.cogs.indica',
    'bot.cogs.sativa',
]
DISCORD_COMMAND_PREFIX = '!'
DISCORD_INTENTS = 'discord.Intents.all()'
DISCORD_MODERATION_WARNING = 'You have been warned.'
DISCORD_OWNER_ID = 154749533429956608
DISCORD_ROLE_PASS = True
DISCORD_TESTING_GUILD_ID = 1300517536001036348

LOGGING_LEVEL = 'INFO'

# OpenAI Chat
OPENAI_CHAT_ADD_COMPLETION_TO_HISTORY = True
OPENAI_CHAT_HEADERS = {
    'Content-Type': 'application/json',
    'OpenAI-Organization': 'org-3LYwtg7DSFJ7RLn9bfk4hATf',
    'User-Agent': 'brandongrahamcobb@icloud.com',
    'OpenAI-Project': 'proj_u5htBCWX0LSHxkw45po1Vfz9',
}
OPENAI_CHAT_MODEL_OUTPUT_LIMITS = {
    'gpt-3.5-turbo': 4096,
    'gpt-4': 8192,
    'gpt-4-32k': 32768,
    'gpt-4o': 4096,         # Initially capped at 4,096; updated to 16,384 in later versions
    'gpt-4o-mini': 16384,
    'gpt-4-turbo': 4096,
    'o1-preview': 32768,
    'o1-mini': 65536,
}
OPENAI_MODEL_CONTEXT_LIMITS = {
    'gpt-3.5-turbo': 4096,
    'gpt-4': 8192,
    'gpt-4-32k': 32768,
    'gpt-4o': 128000,
    'gpt-4o-mini': 128000,
    'gpt-4-turbo': 128000,
    'o1-preview': 128000,
    'o1-mini': 128000,
}
OPENAI_CHAT_MODELS = {
    'current': ['chatgpt-4o-mini-latest', 'o1-preview', 'o1-mini'],
    'deprecated': ['gpt-3.5-turbo', 'gpt-4', 'gpt-4-32k', 'gpt-4o', 'gpt-4o-mini', 'gpt-4-turbo', 'chatgpt-4o-latest'],
}

#OpenAI Moderations
OPENAI_CHAT_MAX_TOKENS = 2000
OPENAI_CHAT_MODERATION = True
OPENAI_CHAT_MODERATION_N = 1
OPENAI_CHAT_MODERATION_MAX_TOKENS = 2000
OPENAI_CHAT_MODERATION_MODEL = 'gpt-4o-mini'
OPENAI_CHAT_MODERATION_RESPONSE_FORMAT = {
'type': 'json_schema',
'json_schema': {
    "name": "moderation",
    "description": "A function that returns moderation results according to a specified schema.",
    "schema": {
      "type": "object",
      "properties": {
        "id": {"type": "string"},
        "model": {"type": "string"},
        "results": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "flagged": {"type": "boolean"},
              "categories": {
                "type": "object",
                "properties": {
                  "sexual": {"type": "boolean"},
                  "hate": {"type": "boolean"},
                  "harassment": {"type": "boolean"},
                  "self-harm": {"type": "boolean"},
                  "sexual/minors": {"type": "boolean"},
                  "hate/threatening": {"type": "boolean"},
                  "violence/graphic": {"type": "boolean"},
                  "self-harm/intent": {"type": "boolean"},
                  "self-harm/instructions": {"type": "boolean"},
                  "harassment/threatening": {"type": "boolean"},
                  "violence": {"type": "boolean"},
                  "carnism": {"type": "boolean"}
                },
                "required": [
                  "sexual",
                  "hate",
                  "harassment",
                  "self-harm",
                  "sexual/minors",
                  "hate/threatening",
                  "violence/graphic",
                  "self-harm/intent",
                  "self-harm/instructions",
                  "harassment/threatening",
                  "violence",
                  "carnism"
                ]
              },
              "category_scores": {
                "type": "object",
                "properties": {
                  "sexual": {"type": "number"},
                  "hate": {"type": "number"},
                  "harassment": {"type": "number"},
                  "self-harm": {"type": "number"},
                  "sexual/minors": {"type": "number"},
                  "hate/threatening": {"type": "number"},
                  "violence/graphic": {"type": "number"},
                  "self-harm/intent": {"type": "number"},
                  "self-harm/instructions": {"type": "number"},
                  "harassment/threatening": {"type": "number"},
                  "violence": {"type": "number"},
                  "carnism": {"type": "number"}
                },
                "required": [
                  "sexual",
                  "hate",
                  "harassment",
                  "self-harm",
                  "sexual/minors",
                  "hate/threatening",
                  "violence/graphic",
                  "self-harm/intent",
                  "self-harm/instructions",
                  "harassment/threatening",
                  "violence",
                  "carnism"
                ]
              }
            },
            "required": ["flagged", "categories", "category_scores"]
          }
        }
      },
      'additionalProperties': False,
      "required": ["id", "model", "results"]
    }
  }
}
OPENAI_CHAT_MODERATION_STOP = ''
OPENAI_CHAT_MODERATION_STORE = False
OPENAI_CHAT_MODERATION_STREAM = False
OPENAI_CHAT_MODERATION_SYS_INPUT = 'You are a JSON moderation assistant.'
OPENAI_CHAT_MODERATION_TEMPERATURE = 1.0
OPENAI_CHAT_MODERATION_TOP_P = 1.0
OPENAI_CHAT_MODEL = 'gpt-4o-mini'
OPENAI_CHAT_N = 1
OPENAI_CHAT_RESPONSE_FORMAT = {
'type': 'json_schema',
'json_schema': {
    "name": "colorize",
    "description": "A function that returns color values for a given request.",
    "schema": {
      "type": "object",
      "properties": {
        "r": {
          "type": "integer",
          "minimum": 0,
          "maximum": 255
        },
        "g": {
          "type": "integer",
          "minimum": 0,
          "maximum": 255
        },
        "b": {
          "type": "integer",
          "minimum": 0,
          "maximum": 255
        }
      },
      "required": ["r", "g", "b"],
      "additionalProperties": False
    }
  }
}
OPENAI_CHAT_MODERATION_USE_HISTORY = False
OPENAI_CHAT_MODERATION_ADD_COMPLETION_TO_HISTORY = False
OPENAI_CHAT_STOP = ''
OPENAI_CHAT_STORE = False
OPENAI_CHAT_STREAM = False
OPENAI_CHAT_SYS_INPUT = ''
OPENAI_CHAT_TOP_P = 1
OPENAI_CHAT_TEMPERATURE = 0.7
OPENAI_CHAT_USE_HISTORY = True
OPENAI_CHAT_USER = 'Brandon Graham Cobb'

OPENAI_FINE_TUNING_RESPONSE_FORMAT = {
  "type": "json_schema",
  "json_schema": {
    "name": "animal_rights_identification",
    "description": "A schema to identify if content qualifies as animal rights activism",
    "schema": {
      "type": "object",
      "properties": {
        "id": {
          "type": "string",
          "description": "Unique identifier for the content being evaluated"
        },
        "results": {
          "type": "boolean",
          "description": "Indicates whether the content qualifies as animal rights activism",
          "enum": [True]
        }
      },
      "required": ["id", "results"]
    }
  }
}

# OpenAI
OPENAI_ENDPOINT_URLS = {
    'audio': 'https://api.openai.com/v1/audio/speech',
    'batch': 'https://api.openai.com/v1/audio/batches',
    'chat': 'https://api.openai.com/v1/chat/completions',
    'embeddings': 'https://api.openai.com/v1/embeddings',
    'files': 'https://api.openai.com/v1/files',
    'fine-tuning': 'https://api.openai.com/v1/fine_tuning/jobs',
    'images': 'https://api.openai.com/v1/images/generations',
    'models': 'https://api.openai.com/v1/models',
    'moderations': 'https://api.openai.com/v1/moderations',
    'uploads': 'https://api.openai.com/v1/uploads',
}

OPENAI_MODERATION_MODEL = 'omni-moderation-latest'
OPENAI_MODERATION_IMAGE = True

OPENAI_TRAINING_FILE = join(PATH_HOME, 'Downloads', 'training_temp.jsonl')

# Scripture Headers
SCRIPTURE_HEADERS = {
    'User-Agent': 'brandongrahamcobb@icloud.com',
    'api-key': '2eb327f99245cd3d68da55370656d6e2'
}
# MySQL
DATABASE_URL = ''

USER_AGENT = 'https://github.com/brandongrahamcobb/py.git'
VERSION = '1.0.0'
