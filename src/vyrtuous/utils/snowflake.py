''' snowflake.py The purpose of this program is to provide the Snowflake utility class.

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
from discord.ext import commands
from discord import app_commands
import discord
import re

class Snowflake:
    def __init__(self, snowflake: int):
        self.snowflake = snowflake
    
    @property
    def snowflake(self) -> int:
        return self._snowflake
    
    @snowflake.setter
    def snowflake(self, snowflake: int):
        if not re.search(r'\d{12,19}', str(snowflake)):
            raise commands.BadArgument(f'The field `{type(self).__name__.lower()}` is invalid.')
        self._snowflake = snowflake

class Converter(commands.Converter):
    
    def __init__(self, snowflake=Snowflake):
        self.snowflake = snowflake
    
    async def convert(self, ctx: commands.Context, arg):
        return self.snowflake(arg).snowflake

class Transformer(app_commands.Transformer):

    def __init__(self, snowflake_cls=Snowflake):
        self.snowflake_cls = snowflake_cls

    async def transform(self, interaction: discord.Interaction, arg):
        return self.snowflake_cls(arg).snowflake

class ChannelSnowflake(Converter):
    def __init__(self):
        super().__init__(Channel)

class GuildSnowflake(Converter):
    def __init__(self):
        super().__init__(Guild)

class MemberSnowflake(Converter):
    def __init__(self):
        super().__init__(Member)

class MessageSnowflake(Converter):
    def __init__(self):
        super().__init__(Message)

class RoleSnowflake(Converter):
    def __init__(self):
        super().__init__(Role)

class AppChannelSnowflake(Transformer):
    def __init__(self):
        super().__init__(Channel)

class AppGuildSnowflake(Transformer):
    def __init__(self):
        super().__init__(Guild)

class AppMemberSnowflake(Transformer):
    def __init__(self):
        super().__init__(Member)

class AppMessageSnowflake(Transformer):
    def __init__(self):
        super().__init__(Message)

class AppRoleSnowflake(Transformer):
    def __init__(self):
        super().__init__(Role)

class Channel(Snowflake):
    pass

class Guild(Snowflake):
    pass

class Member(Snowflake):
    pass

class Message(Snowflake):
    pass

class Role(Snowflake):
    pass