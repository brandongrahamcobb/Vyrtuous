''' config.py  The purpose of this program is to provide my primary configuaration script from cd ../.
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
from os import makedirs
from os.path import isfile, dirname
from typing import Dict, Any
from .load_yaml import load_yaml
from .setup_logging import logger
from .prompt_for_values import prompt_for_values
from .helpers import *
import yaml

class Config:
    _config = None  # Class variable to hold the config config

    @classmethod
    def get_config(cls) -> Dict[str, Any]:
        logger.info('Attempting to load configuration.')
        
        if cls._config is None:
            if isfile(PATH_CONFIG_YAML):
                logger.info(f'Config file found at {PATH_CONFIG_YAML}. Loading configuration.')
                config = load_yaml(PATH_CONFIG_YAML)
                
                # Ensure 'api_keys' exists and set defaults
                config['api_keys'] = config.get('api_keys', {})
                logger.info('API keys initialized.')

                # Loop through and prompt for each API key and associated values
                for i in range(1, 21):
                    key = f'api_key_{i}'
                    logger.info(f'Prompting for values for {key}.')
                    config['api_keys'][key] = config['api_keys'].get(key, {})
                    config['api_keys'][key]['api_key'] = prompt_for_values(f'Enter API key {i}', config['api_keys'][key].get('api_key', ''))
                    config['api_keys'][key]['client_id'] = prompt_for_values(f'Enter client ID for API key {i}', config['api_keys'][key].get('client_id', ''))
                    config['api_keys'][key]['client_secret'] = prompt_for_values(f'Enter client secret for API key {i}', config['api_keys'][key].get('client_secret', ''))
                    config['api_keys'][key]['redirect_uri'] = prompt_for_values(f'Enter redirect URI for API key {i}', config['api_keys'][key].get('redirect_uri', ''))

                # Prompt for additional configuration values
                logger.info('Prompting for other configuration values.')
                config['discord_character_limit'] = prompt_for_values('Discord character limit?', config.get('discord_character_limit', DISCORD_CHARACTER_LIMIT))
                config['discord_command_prefix'] = prompt_for_values('Discord command prefix?', config.get('discord_command_prefix', DISCORD_COMMAND_PREFIX))
                config['discord_moderation_warning'] = prompt_for_values('What should be sent to users if their message was moderated?', config.get('discord_moderation_warning', DISCORD_MODERATION_WARNING))
                config['discord_owner_id'] = prompt_for_values('Discord Owner ID?', config.get('discord_owner_id', DISCORD_OWNER_ID))
                config['discord_role_pass'] = prompt_for_values('Should vegans not be moderated?', config.get('discord_role_pass', DISCORD_ROLE_PASS))
                config['discord_testing_guild_id'] = prompt_for_values('What is the Discord testing guild ID?', config.get('discord_testing_guild_id', DISCORD_TESTING_GUILD_ID))
                config['discord_token'] = prompt_for_values('What is the Discord token?', config.get('discord_token', ''))
                config['logging_level'] = prompt_for_values('What is the logging level (DEBUG, INFO, etc.)?', config.get('logging_level', LOGGING_LEVEL))
                config['openai_chat_add_completion_to_history'] = prompt_for_values('Should completions be added to conversations?', config.get('openai_chat_add_completion_to_history', OPENAI_CHAT_ADD_COMPLETION_TO_HISTORY))
                config['openai_chat_max_tokens'] = prompt_for_values('What is the max tokens per OpenAI chat completion?', config.get('openai_chat_max_tokens', OPENAI_CHAT_MAX_TOKENS))
                config['openai_chat_model'] = prompt_for_values("Which chat model would you like to use for OpenAI's ChatGPT?", config.get('openai_chat_model', OPENAI_CHAT_MODEL))
                config['openai_chat_moderation_model'] = prompt_for_values('Which OpenAI completions model would you like to use for moderation?', config.get('openai_chat_moderation_model', OPENAI_CHAT_MODERATION_MODEL))
                config['openai_chat_moderation'] = prompt_for_values('Enable or disable OpenAI text moderation (True/False)?', config.get('openai_moderation_text', OPENAI_CHAT_MODERATION))
                config['openai_chat_store'] = prompt_for_values('Store OpenAI completions (True/False)?', config.get('openai_chat_store', OPENAI_CHAT_STORE))
                config['openai_chat_stream'] = prompt_for_values('Enable or disable OpenAI completions streaming (True/False)?', config.get('openai_chat_stream', OPENAI_CHAT_STREAM))
                config['openai_chat_stop'] = prompt_for_values('What might be the OpenAI stop criteria for completions?', config.get('openai_chat_stop', OPENAI_CHAT_STOP))
                default_chat_sys_input = OPENAI_CHAT_SYS_INPUT if OPENAI_CHAT_MODEL in OPENAI_CHAT_MODELS['deprecated'] else ''
                config['openai_chat_sys_input'] = prompt_for_values('What is the OpenAI completions system input?', config.get('openai_chat_sys_input', default_chat_sys_input))
                config['openai_chat_temperature'] = prompt_for_values('What is the OpenAI completions temperature (0.0 to 2.0)?', config.get('openai_chat_temperature', OPENAI_CHAT_TEMPERATURE))
                config['openai_chat_top_p'] = prompt_for_values('What should the top p be for OpenAI completions?', config.get('openai_chat_top_p', OPENAI_CHAT_TOP_P))
                config['openai_chat_use_history'] = prompt_for_values('Should OpenAI moderations use history?', config.get('openai_chat_use_history', OPENAI_CHAT_USE_HISTORY))
                config['openai_chat_user'] = prompt_for_values('What is your OpenAI username?', config.get('openai_chat_user', OPENAI_CHAT_USER))
                config['openai_moderation_image'] = prompt_for_values('Enable or disable OpenAI image moderation (True/False)?', config.get('openai_moderation_image', OPENAI_MODERATION_IMAGE))
                config['openai_moderation_model'] = prompt_for_values('Which model do you want for OpenAI image moderation?', config.get('openai_moderation_model', OPENAI_MODERATION_MODEL))
                config['openai_organization'] = prompt_for_values('What is the OpenAI-Organization ID?', config.get('openai_organization', OPENAI_CHAT_HEADERS['OpenAI-Organization']))
                config['openai_project'] = prompt_for_values('What is the OpenAI-Project ID?', config.get('openai_project', OPENAI_CHAT_HEADERS['OpenAI-Project']))
                config['user_agent'] = prompt_for_values('What should be the User-Agent header?', config.get('user_agent', USER_AGENT))
                config['version'] = prompt_for_values('Would you like to override the bot version?', config.get('version', VERSION))

                logger.info('Configuration loaded and all values prompted.')
            else:
                logger.info(f'Config file not found. Creating new config file at {PATH_CONFIG_YAML}.')
                makedirs(dirname(PATH_CONFIG_YAML), exist_ok=True)
                config = {
                    'api_keys': {
                        f'api_key_{i}': {
                            'api_key': prompt_for_values(f'Enter API key {i}', ''),
                            'client_id': prompt_for_values(f'Enter client ID for API key {i}', ''),
                            'client_secret': prompt_for_values(f'Enter client secret for API key {i}', ''),
                            'redirect_uri': prompt_for_values(f'Enter redirect URI for API key {i}', '')
                        } for i in range(1, 21)
                    },
                    'discord_character_limit': prompt_for_values('Discord character limit?', DISCORD_CHARACTER_LIMIT),
                    'discord_command_prefix': prompt_for_values('Discord command prefix?', DISCORD_COMMAND_PREFIX),
                    'discord_moderation_warning': prompt_for_values('Message for moderated users?', DISCORD_MODERATION_WARNING),
                    'discord_owner_id': prompt_for_values('Discord Owner ID?', DISCORD_OWNER_ID),
                    'discord_role_pass': prompt_for_values('Should vegans not be moderated?', DISCORD_ROLE_PASS),
                    'discord_testing_guild_id': prompt_for_values('Discord testing guild ID?', DISCORD_TESTING_GUILD_ID),
                    'discord_token': prompt_for_values('Discord token?', ''),
                    'logging_level': prompt_for_values('Logging level (DEBUG, INFO, etc.)?', LOGGING_LEVEL),
                    'openai_chat_add_completion_to_history': prompt_for_values('Add completions to conversations?', OPENAI_CHAT_ADD_COMPLETION_TO_HISTORY),
                    'openai_chat_max_tokens': prompt_for_values('Max tokens per OpenAI chat completion?', OPENAI_CHAT_MAX_TOKENS),
                    'openai_chat_moderation_model': prompt_for_values('Chat model for moderation?', OPENAI_CHAT_MODERATION_MODEL),
                    'openai_chat_model': prompt_for_values('Chat model for OpenAI ChatGPT?', OPENAI_CHAT_MODEL),
                    'openai_chat_moderation': prompt_for_values('Enable or disable OpenAI text moderation (True/False)?', OPENAI_CHAT_MODERATION),
                    'openai_chat_store': prompt_for_values('Store OpenAI completions (True/False)?', OPENAI_CHAT_STORE),
                    'openai_chat_stream': prompt_for_values('Enable or disable OpenAI completions streaming (True/False)?', OPENAI_CHAT_STREAM),
                    'openai_chat_stop': prompt_for_values('OpenAI stop criteria for completions?', OPENAI_CHAT_STOP),
                    'openai_chat_sys_input': prompt_for_values('OpenAI completions system input?', OPENAI_CHAT_SYS_INPUT),
                    'openai_chat_temperature': prompt_for_values('OpenAI completions temperature (0.0 to 2.0)?', OPENAI_CHAT_TEMPERATURE),
                    'openai_chat_top_p': prompt_for_values('Top p for OpenAI completions?', OPENAI_CHAT_TOP_P),
                    'openai_chat_use_history': prompt_for_values('Use chat history?', OPENAI_CHAT_USE_HISTORY),
                    'openai_chat_user': prompt_for_values('Your OpenAI username?', OPENAI_CHAT_USER),
                    'openai_moderation_image': prompt_for_values('Enable or disable OpenAI image moderation (True/False)?', OPENAI_MODERATION_IMAGE),
                    'openai_moderation_model': prompt_for_values('Model for OpenAI image moderation?', OPENAI_MODERATION_MODEL),
                    'openai_organization': prompt_for_values('OpenAI-Organization ID?', OPENAI_CHAT_HEADERS['OpenAI-Organization']),
                    'openai_project': prompt_for_values('OpenAI-Project ID?', OPENAI_CHAT_HEADERS['OpenAI-Project']),
                    'user_agent': prompt_for_values('User-Agent header?', USER_AGENT),
                    'version': prompt_for_values('Bot version override?', VERSION),
                }

                with open(PATH_CONFIG_YAML, 'w') as file:
                    yaml.dump(config, file)
                    logger.info(f'Configuration saved to {PATH_CONFIG_YAML}.')

            cls._config = config
            logger.info('Configuration successfully loaded.')
        else:
            logger.info('Configuration already loaded, returning cached version.')
        
        return cls._config
