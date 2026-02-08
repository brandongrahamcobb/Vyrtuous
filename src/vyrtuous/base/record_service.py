"""!/bin/python3

service.py The purpose of this program is to the parent class to all service classes.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from typing import Union

from discord.ext import commands
import discord

from vyrtuous.commands.messaging.state_service import StateService


class RecordService:
    @classmethod
    async def enforce_or_undo(
        cls,
        ctx,
        source: Union[commands.Context, discord.Interaction, discord.Message],
        state: StateService,
    ):
        obj = await ctx.record.select(**ctx.source_kwargs, singular=True)
        if obj:
            await cls.undo(ctx=ctx, source=source, state=state)
        else:
            await cls.enforce(ctx=ctx, source=source, state=state)
