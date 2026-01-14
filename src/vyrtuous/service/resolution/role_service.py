"""channel_service.py The purpose of this program is to provide the channel_service module.
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
"""

from typing import Optional, Union

import discord

from vyrtuous.service.logging_service import logger


async def resolve_role(
    ctx_interaction_or_message,
    role_str: Optional[Union[int, str, discord.Role]],
) -> discord.Role:
    try:
        if isinstance(role_str, discord.Role):
            logger.info(f"Direct role: {role_str.id}")
            return role_str
        if isinstance(role_str, int):
            role = ctx_interaction_or_message.guild.get_role(role_str)
            if role:
                logger.info(f"Resolved role by int ID: {role.id}")
                return role
        if isinstance(role_str, str):
            if role_str.isdigit():
                role = ctx_interaction_or_message.guild.get_role(int(role_str))
                if role:
                    logger.info(f"Resolved role by str ID: {role.id}")
                    return role
            if role_str.startswith("<@&") and role_str.endswith(">"):
                role_id = int(role_str[3:-1])
                role = ctx_interaction_or_message.guild.get_role(role_id)
                if role:
                    logger.info(f"Role mention resolved: {role.id}")
                    return role
    except Exception as e:
        logger.warning(f"Role resolution error: {str(e).capitalize()}")
        raise
    raise ValueError("Role could not be resolved from role_str `{role_str}`.")
