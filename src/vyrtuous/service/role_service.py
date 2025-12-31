''' channel_service.py The purpose of this program is to provide the channel_service module.
    Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

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
from typing import Optional, Union
from vyrtuous.utils.setup_logging import logger
import discord

class RoleService:
    
    async def resolve_role(
        self,
        ctx_interaction_or_message,
        scope: Optional[Union[int, str, discord.Role]]
    ) -> discord.Role:
        try:
            if isinstance(scope, discord.Role):
                logger.debug(f'Direct role: {scope.id}')
                return scope
            if isinstance(scope, int):
                role = ctx_interaction_or_message.guild.get_role(scope)
                if role:
                    logger.debug(f'Resolved role by int ID: {role.id}')
                    return role
            if isinstance(scope, str):
                if scope.isdigit():
                    role = ctx_interaction_or_message.guild.get_role(int(scope))
                    if role:
                        logger.debug(f'Resolved role by str ID: {role.id}')
                        return role
                if scope.startswith('<@&') and scope.endswith('>'):
                    role_id = int(scope[3:-1])
                    role = ctx_interaction_or_message.guild.get_role(role_id)
                    if role:
                        logger.debug(f'Role mention resolved: {role.id}')
                        return role
        except Exception as e:
            logger.warning(f'Role resolution error: {str(e).capitalize()}')
            raise
        raise ValueError('Role could not be resolved from scope `{scope}`.')
