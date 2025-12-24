from discord.ext import commands
from discord import app_commands
import discord
import re

class Snowflake:
    def __init__(self, snowflake: int):
        self.snowflake = snowflake

    def __int__(self):
        return self._snowflake
    
    @property
    def snowflake(self) -> int:
        return self._snowflake
    
    @snowflake.setter
    def snowflake(self, snowflake: int):
        if not re.search(r"\d{12,19}", str(snowflake)):
            raise commands.BadArgument(f"{type(self).__name__} is invalid")
        self._snowflake = snowflake

class Converter(commands.Converter):
    
    def __init__(self, snowflake=Snowflake):
        self.snowflake = snowflake
    
    async def convert(self, ctx: commands.Context, arg):
        return self.snowflake(arg)

class Transformer(app_commands.Transformer):

    def __init__(self, snowflake_cls=Snowflake):
        self.snowflake_cls = snowflake_cls

    async def transform(self, interaction: discord.Interaction, arg):
        return self.snowflake(arg)

class ChannelSnowflake(Converter):
    def __init__(self):
        super().__init__(Snowflake)

class GuildSnowflake(Converter):
    def __init__(self):
        super().__init__(Snowflake)

class MemberSnowflake(Converter):
    def __init__(self):
        super().__init__(Snowflake)

class MessageSnowflake(Converter):
    def __init__(self):
        super().__init__(Snowflake)

class RoleSnowflake(Converter):
    def __init__(self):
        super().__init__(Snowflake)

class AppChannelSnowflake(Transformer):
    def __init__(self):
        super().__init__(Snowflake)

class AppGuildSnowflake(Transformer):
    def __init__(self):
        super().__init__(Snowflake)

class AppMemberSnowflake(Transformer):
    def __init__(self):
        super().__init__(Snowflake)

class AppMessageSnowflake(Transformer):
    def __init__(self):
        super().__init__(Snowflake)

class AppRoleSnowflake(Transformer):
    def __init__(self):
        super().__init__(Snowflake)