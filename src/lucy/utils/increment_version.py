''' increment_version.py  The purpose of this program is to provide persistent versioning from cd ../.
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
from typing import Any, Dict
from utils.setup_logging import logger

import yaml

def increment_version(config: Dict[str, Any], path_config_yaml):
    try:
        logger.info('Starting version increment process.')

        # Retrieve the current version
        current_version = config.get('version', '0.0.0')
        logger.debug(f'Current version: {current_version}')

        # Parse and increment the version
        major, minor, patch = map(int, current_version.split('.'))
        patch += 1
        if patch >= 10:
            patch = 0
            minor += 1
        if minor >= 10:
            minor = 0
            major += 1

        new_version = f'{major}.{minor}.{patch}'
        logger.info(f'New version generated: {new_version}')

        # Update the version in the config
        config['version'] = new_version

        # Write the updated config back to the YAML file
        with open(path_config_yaml, 'w') as file:
            yaml.dump(config, file)
        logger.info(f'Version updated successfully in the YAML file: {path_config_yaml}')

    except Exception as e:
        logger.error(f'An error occurred during version increment: {e}')
        raise
