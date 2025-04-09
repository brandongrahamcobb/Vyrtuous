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
from lucy.utils.load_contents import load_contents
from lucy.utils.setup_logging import logger
from os.path import dirname, abspath, expanduser, join

import discord

# Base and Home Paths
DIR_BASE = abspath(join(dirname(dirname(dirname(__file__)))))
DIR_HOME = expanduser('~')
DIR_TEMP = join(DIR_BASE, 'lucy', 'temp')

# Script Paths
PATH_ARPP = join(DIR_BASE, 'lucy', 'utils', 'api_request_parallel_processor.py')
PATH_AVERAGE_SCORE = join(DIR_BASE, 'lucy', 'utils', 'average_score.py')
PATH_BENCHMARK = join(DIR_BASE, 'lucy', 'utils', 'benchmark.py')
PATH_CHEMISTRY = join(DIR_BASE, 'lucy', 'utils', 'chemistry.py')
PATH_CITATION_MANAGER = join(DIR_BASE, 'lucy', 'utils', 'citation_manager.py')
PATH_CLEAR_SCREEN = join(DIR_BASE, 'lucy', 'utils', 'clear_screen.py')
PATH_CONFIG = join(DIR_BASE, 'lucy', 'utils', 'config.py')
PATH_CONFIG_YAML = join(DIR_BASE, 'lucy', '.config', 'config.yaml')
PATH_DISCORD_BOT = join(DIR_BASE, 'lucy', 'bots', 'discord_bot.py')
PATH_DISCORD_OAUTH = join(DIR_BASE, 'lucy', 'utils', 'discord_oauth.py')
PATH_FINE_TUNING = join(DIR_BASE, 'lucy', 'utils', 'fine_tuning.py')
PATH_FONT = join(DIR_HOME, 'Vyrtuous', 'Roboto-Regular.ttf')
PATH_FORMAT_ERROR_CHECK = join(DIR_BASE, 'lucy', 'utils', 'format_error_check.py')
PATH_FRAMES = join(DIR_BASE, 'lucy', 'utils', 'frames.py')
PATH_GAME = join(DIR_BASE, 'lucy', 'utils', 'game.py')
PATH_GET_SCRIPTURE = join(DIR_BASE, 'lucy', 'utils', 'get_scripture.py')
PATH_GOOGLE = join(DIR_BASE, 'lucy', 'utils', 'google.py')
PATH_HELPERS = join(DIR_BASE, 'lucy', 'utils', 'helpers.py')
PATH_HYBRID = join(DIR_BASE, 'lucy', 'cogs', 'hybrid.py')
PATH_IMAGE = join(DIR_BASE, 'lucy', 'utils', 'image.py')
PATH_INCREMENT_VERSION = join(DIR_BASE, 'lucy', 'utils', 'increment_version.py')
PATH_INDICA = join(DIR_BASE, 'lucy', 'cogs', 'indica.py')
PATH_LINKEDIN_BOT = join(DIR_BASE, 'lucy', 'bots', 'linkedin_bot.py')
PATH_LINKEDIN_OAUTH = join(DIR_BASE, 'lucy', 'utils', 'linkedin_oauth.py')
PATH_LOAD_CONTENTS = join(DIR_BASE, 'lucy', 'utils', 'load_contents.py')
PATH_LOAD_YAML = join(DIR_BASE, 'lucy', 'utils', 'load_yaml.py')
PATH_LOG = join(DIR_BASE, 'lucy', '.log', 'discord.log')
PATH_MAIN = join(DIR_BASE, 'lucy', 'main.py')
PATH_MESSAGE = join(DIR_BASE, 'lucy', 'utils', 'message.py')
PATH_NLP_UTILS = join(DIR_BASE, 'lucy', 'utils', 'nlp_utils.py')
PATH_OPENAI_REQUESTS = join(DIR_HOME, 'Downloads', 'queued_requests.json')
PATH_OPENAI_RESULTS = join(DIR_HOME, 'Downloads', 'processed_results.json')
PATH_PAGINATOR = join(DIR_BASE, 'lucy', 'utils', 'paginator.py')
PATH_PDF_MANAGER = join(DIR_BASE, 'lucy', 'utils', 'pdf_manager.py')
PATH_PATREON_OAUTH = join(DIR_BASE, 'lucy', 'utils', 'patreon_oauth.py')
PATH_PREDICATOR = join(DIR_BASE, 'lucy', 'utils', 'predicator.py')
PATH_PROMPT_FOR_VALUES = join(DIR_BASE, 'lucy', 'utils', 'prompt_for_values.py')
PATH_ROLE_MANAGER = join(DIR_BASE, 'lucy', 'utils', 'role_manager.py')
PATH_SATIVA = join(DIR_BASE, 'lucy', 'cogs', 'sativa.py')
PATH_SCRIPT = join(DIR_BASE, 'lucy', 'utils', 'script.py')
PATH_SETUP_LOGGING = join(DIR_BASE, 'lucy', 'utils', 'setup_logging.py')
PATH_TAG = join(DIR_BASE, 'lucy', 'utils', 'tag.py')
PATH_TOML = join(DIR_HOME, 'Vyrtuous', 'pyproject.toml')
PATH_TRAINING = join(DIR_HOME, '..', 'training.jsonl')
PATH_TWITCH_BOT = join(DIR_BASE, 'lucy', 'bots', 'twitch_bot.py')
PATH_TWITCH_OAUTH = join(DIR_BASE, 'lucy', 'utils', 'twitch_oauth.py')
PATH_UNIQUE_PAIRS = join(DIR_BASE, 'lucy', 'utils', 'unique_pairs.py')
PATH_USERS = join(DIR_BASE, 'lucy', '.users', 'users.yaml')

# Main
CONTENTS_MAIN = load_contents(PATH_MAIN)

# Utilities
CONTENTS_ARPP = load_contents(PATH_ARPP)
CONTENTS_AVERAGE_SCORE = load_contents(PATH_AVERAGE_SCORE)
CONTENTS_BENCHMARK = load_contents(PATH_BENCHMARK)
CONTENTS_CLEAR_SCREEN = load_contents(PATH_CLEAR_SCREEN)
CONTENTS_CHEMISTRY = load_contents(PATH_CHEMISTRY)
CONTENTS_DISCORD_BOT = load_contents(PATH_DISCORD_BOT)
CONTENTS_DISCORD_OAUTH = load_contents(PATH_DISCORD_OAUTH)
CONTENTS_FINE_TUNING = load_contents(PATH_FINE_TUNING)
CONTENTS_FORMAT_ERROR_CHECK = load_contents(PATH_FORMAT_ERROR_CHECK)
CONTENTS_FRAMES = load_contents(PATH_FRAMES)
CONTENTS_GAME = load_contents(PATH_GAME)
CONTENTS_GOOGLE = load_contents(PATH_GOOGLE)
CONTENTS_IMAGE = load_contents(PATH_IMAGE)
CONTENTS_HELPERS = load_contents(PATH_HELPERS)
CONTENTS_INCREMENT_VERSION = load_contents(PATH_INCREMENT_VERSION)
CONTENTS_LINKEDIN_BOT = load_contents(PATH_LINKEDIN_BOT)
CONTENTS_LINKEDIN_OAUTH = load_contents(PATH_LINKEDIN_OAUTH)
CONTENTS_LOAD_CONTENTS= load_contents(PATH_LOAD_CONTENTS)
CONTENTS_LOAD_YAML = load_contents(PATH_LOAD_YAML)
CONTENTS_MESSAGE = load_contents(PATH_MESSAGE)
CONTENTS_PROMPT_FOR_VALUES = load_contents(PATH_PROMPT_FOR_VALUES)
CONTENTS_PAGINATOR = load_contents(PATH_PAGINATOR)
CONTENTS_PATREON_OAUTH = load_contents(PATH_PATREON_OAUTH)
CONTENTS_PDF_MANAGER = load_contents(PATH_PDF_MANAGER)
CONTENTS_PREDICATOR = load_contents(PATH_PREDICATOR)
CONTENTS_SCRIPT = load_contents(PATH_SCRIPT)
CONTENTS_SETUP_LOGGING = load_contents(PATH_SETUP_LOGGING)
CONTENTS_TAG = load_contents(PATH_TAG)
CONTENTS_TWITCH_BOT = load_contents(PATH_TWITCH_BOT)
CONTENTS_TWITCH_OAUTH = load_contents(PATH_TWITCH_OAUTH)
CONTENTS_UNIQUE_PAIRS = load_contents(PATH_UNIQUE_PAIRS)
SUM_OF_UTILITIES = f'''
    {CONTENTS_ARPP} and {CONTENTS_AVERAGE_SCORE} and {CONTENTS_BENCHMARK} and {CONTENTS_CLEAR_SCREEN} \n
    {CONTENTS_DISCORD_BOT} and {CONTENTS_DISCORD_OAUTH} \n
    {CONTENTS_FINE_TUNING} and {CONTENTS_FORMAT_ERROR_CHECK} and {CONTENTS_FRAMES} \n
    {CONTENTS_GAME} and {CONTENTS_GOOGLE} and {CONTENTS_HELPERS} and {CONTENTS_INCREMENT_VERSION} \n
    {CONTENTS_LINKEDIN_BOT} and {CONTENTS_LINKEDIN_OAUTH} and {CONTENTS_LOAD_CONTENTS} and {CONTENTS_LOAD_YAML} and {CONTENTS_MESSAGE} and {CONTENTS_PAGINATOR} and {CONTENTS_PATREON_OAUTH} and {CONTENTS_PDF_MANAGER} AND {CONTENTS_PREDICATOR} and {CONTENTS_PROMPT_FOR_VALUES}
    {CONTENTS_SCRIPT} and {CONTENTS_SETUP_LOGGING} and {CONTENTS_TAG} and {CONTENTS_TWITCH_BOT} {CONTENTS_TWITCH_OAUTH} {CONTENTS_UNIQUE_PAIRS} \n
'''
# MySQL
DATABASE_URL = ''

# Discord
CONTENTS_HYBRID = load_contents(PATH_HYBRID)
CONTENTS_INDICA = load_contents(PATH_INDICA)
CONTENTS_SATIVA = load_contents(PATH_SATIVA)
SUM_OF_COGS = f'''
    {CONTENTS_HYBRID} and {CONTENTS_INDICA} and {CONTENTS_SATIVA}
'''
DISCORD_CHARACTER_LIMITS = [2000, 4000]
DISCORD_CHARACTER_LIMIT = 2000
DISCORD_COGS = [
    'lucy.cogs.hybrid',
    'lucy.cogs.indica',
    'lucy.cogs.ruderalis',
    'lucy.cogs.sativa',
]
DISCORD_COMMAND_PREFIX = '!'
DISCORD_INTENTS = discord.Intents.all()
DISCORD_MODERATION_WARNING = 'You have been warned.'
DISCORD_OWNER_ID = 154749533429956608
DISCORD_RELEASE_MODE = False
DISCORD_ROLE_PASS = 1308689505158565918
DISCORD_TESTING_GUILD_ID = 1300517536001036348

FTP_HOSTNAME = ''
FTP_PASSWORD = ''
FTP_PDF_URL = 'http://brandongcobb.com/pdfs'
FTP_PUBLIC_URL = 'http://brandongcobb.com'
FTP_USER = ''


LOGGING_LEVEL = 'INFO'

# OpenAI Chat
OPENAI_CHAT_ADD_COMPLETION_TO_HISTORY = True
OPENAI_CHAT_COLORIZE_RESPONSE_FORMAT = {
'type': 'json_schema',
'json_schema': {
    'name': 'colorize',
    'description': 'A function that returns color values for a given request.',
    'schema': {
      'type': 'object',
      'properties': {
        'r': {
          'type': 'integer',
          'minimum': 0,
          'maximum': 255
        },
        'g': {
          'type': 'integer',
          'minimum': 0,
          'maximum': 255
        },
        'b': {
          'type': 'integer',
          'minimum': 0,
          'maximum': 255
        }
      },
      'required': ['r', 'g', 'b'],
      'additionalProperties': False
    }
  }
}
OPENAI_CHAT_COMPLETION = True
OPENAI_CHAT_HEADERS = {
    'Content-Type': 'application/json',
    'OpenAI-Organization': 'org-3LYwtg7DSFJ7RLn9bfk4hATf',
    'User-Agent': 'brandongrahamcobb@icloud.com',
    'OpenAI-Project': 'proj_u5htBCWX0LSHxkw45po1Vfz9',
}
OPENAI_CHAT_MODELS = {
    'current': ['chatgpt-4o-mini-latest', 'o1-preview', 'o1-mini'],
    'deprecated': ['gpt-3.5-turbo', 'gpt-4', 'gpt-4-32k', 'gpt-4o', 'gpt-4o-mini', 'gpt-4-turbo', 'chatgpt-4o-latest'],
}

#OpenAI Moderations
OPENAI_CHAT_MODERATION = True
OPENAI_CHAT_MODERATION_N = 1
OPENAI_CHAT_MODERATION_MODEL = 'gpt-4o-mini'
OPENAI_CHAT_MODERATION_RESPONSE_FORMAT = {
'type': 'json_schema',
'json_schema': {
    'name': 'moderation',
    'description': 'A function that returns moderation results according to a specified schema.',
    'schema': {
      'type': 'object',
      'properties': {
        'id': {'type': 'string'},
        'model': {'type': 'string'},
        'results': {
          'type': 'array',
          'items': {
            'type': 'object',
            'properties': {
              'flagged': {'type': 'boolean'},
              'categories': {
                'type': 'object',
                'properties': {
                  'sexual': {'type': 'boolean'},
                  'hate': {'type': 'boolean'},
                  'harassment': {'type': 'boolean'},
                  'self-harm': {'type': 'boolean'},
                  'sexual/minors': {'type': 'boolean'},
                  'hate/threatening': {'type': 'boolean'},
                  'violence/graphic': {'type': 'boolean'},
                  'self-harm/intent': {'type': 'boolean'},
                  'self-harm/instructions': {'type': 'boolean'},
                  'harassment/threatening': {'type': 'boolean'},
                  'violence': {'type': 'boolean'},
                  'academic-dishonesty': {'type': 'boolean'},
                  'animal-derived-technology': {'type': 'boolean'}
                },
                'required': [
                  'sexual',
                  'hate',
                  'harassment',
                  'self-harm',
                  'sexual/minors',
                  'hate/threatening',
                  'violence/graphic',
                  'self-harm/intent',
                  'self-harm/instructions',
                  'harassment/threatening',
                  'violence',
                  'academic-dishonesty',
                  'animal-derived-technology'
                ]
              },
              'category_scores': {
                'type': 'object',
                'properties': {
                  'sexual': {'type': 'number'},
                  'hate': {'type': 'number'},
                  'harassment': {'type': 'number'},
                  'self-harm': {'type': 'number'},
                  'sexual/minors': {'type': 'number'},
                  'hate/threatening': {'type': 'number'},
                  'violence/graphic': {'type': 'number'},
                  'self-harm/intent': {'type': 'number'},
                  'self-harm/instructions': {'type': 'number'},
                  'harassment/threatening': {'type': 'number'},
                  'violence': {'type': 'number'},
                  'academic-dishonesty': {'type': 'number'},
                  'animal-derived-technology': {'type': 'boolean'}
                },
                'required': [
                  'sexual',
                  'hate',
                  'harassment',
                  'self-harm',
                  'sexual/minors',
                  'hate/threatening',
                  'violence/graphic',
                  'self-harm/intent',
                  'self-harm/instructions',
                  'harassment/threatening',
                  'violence',
                  'academic-dishonesty',
                  'animal-derived-technology'
                ]
              }
            },
            'required': ['flagged', 'categories', 'category_scores']
          }
        }
      },
      'additionalProperties': False,
      'required': ['id', 'model', 'results']
    }
  }
}
OPENAI_CHAT_MODERATION_STOP = ''
OPENAI_CHAT_MODERATION_STORE = False
OPENAI_CHAT_MODERATION_STREAM = False
OPENAI_CHAT_MODERATION_SYS_INPUT = 'You are a JSON moderation assistant.'
OPENAI_CHAT_MODERATION_TEMPERATURE = 1.0
OPENAI_CHAT_MODERATION_TOP_P = 1.0
OPENAI_CHAT_MODERATION_USE_HISTORY = False
OPENAI_CHAT_MODEL = 'gpt-4o-mini'
OPENAI_CHAT_N = 1
OPENAI_CHAT_RESPONSE_FORMAT = None
OPENAI_CHAT_MODERATION_USE_HISTORY = False
OPENAI_CHAT_MODERATION_ADD_COMPLETION_TO_HISTORY = False
OPENAI_CHAT_STOP = ''
OPENAI_CHAT_STORE = False
OPENAI_CHAT_STREAM = False
OPENAI_CHAT_SYS_INPUT = f'Your main.py is {CONTENTS_MAIN}. Your cogs are {SUM_OF_COGS}. Your utilities {SUM_OF_UTILITIES}. Short conversational responses only.'
OPENAI_CHAT_TOP_P = 1
OPENAI_CHAT_TEMPERATURE = 0.7
OPENAI_CHAT_USE_HISTORY = True
OPENAI_CHAT_USER = 'Brandon Graham Cobb'
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
OPENAI_FINE_TUNING_RESPONSE_FORMAT = {
  'type': 'json_schema',
  'json_schema': {
    'name': 'animal_rights_identification',
    'description': 'A schema to identify if content qualifies as animal rights activism',
    'schema': {
      'type': 'object',
      'properties': {
        'id': {
          'type': 'string',
          'description': 'Unique identifier for the content being evaluated'
        },
        'results': {
          'type': 'boolean',
          'description': 'Indicates whether the content qualifies as animal rights activism',
          'enum': [True]
        }
      },
      'required': ['id', 'results']
    }
  }
}
OPENAI_MODEL_CONTEXT_LIMITS = {
    'ft:gpt-4o-mini-2024-07-18:spawd:vyrtuous:AjZpTNN2':128000,
    'gpt-3.5-turbo': 4096,
    'gpt-4': 8192,
    'gpt-4-32k': 32768,
    'gpt-4o': 128000,
    'gpt-4o-mini': 128000,
    'gpt-4-turbo': 128000,
    'o1-preview': 128000,
    'o1-mini': 128000,
}
OPENAI_MODEL_OUTPUT_LIMITS = {
    'ft:gpt-4o-mini-2024-07-18:spawd:vyrtuous:AjZpTNN2': 16384,
    'gpt-3.5-turbo': 4096,
    'gpt-4': 8192,
    'gpt-4-32k': 32768,
    'gpt-4o': 4096,         # Initially capped at 4,096; updated to 16,384 in later versions
    'gpt-4o-mini': 16384,
    'gpt-4-turbo': 4096,
    'o1-preview': 32768,
    'o1-mini': 16384,
}
OPENAI_MODERATION_MODEL = 'omni-moderation-latest'
OPENAI_MODERATION_IMAGE = True

# Scripture Headers
SCRIPTURE_HEADERS = {
    'User-Agent': 'brandongrahamcobb@icloud.com',
    'api-key': '2eb327f99245cd3d68da55370656d6e2'
}

USER_AGENT = 'https://github.com/brandongrahamcobb/lucy.git'
VERSION = '1.0.0'
