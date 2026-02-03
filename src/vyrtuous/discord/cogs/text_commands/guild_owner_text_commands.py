"""guild_owner_commands.py A discord.py cog containing guild owner commands for the Vyrtuous bot.

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

from discord.ext import commands

from vyrtuous.db.roles.admin.administrator_service import AdministratorRoleService
from vyrtuous.db.roles.owner.guild_owner_service import guild_owner_predicator
from vyrtuous.discord.bot.discord_bot import DiscordBot
from vyrtuous.discord.cogs.help_command import skip_help_discovery
from vyrtuous.discord.discord_object_service import DiscordObject
from vyrtuous.discord.fields.snowflake import MemberSnowflake, RoleSnowflake
from vyrtuous.discord.messaging.message_service import MessageService
from vyrtuous.discord.messaging.state_service import StateService
from vyrtuous.permissions.invincibility import Invincibility


class GuildOwnerTextCommands(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.message_service = MessageService(self.bot)

    @commands.command(name="admin", help="Toggle administrator role.")
    @guild_owner_predicator()
    async def toggle_administrator_by_role_text_command(
        self, ctx: commands.Context, role: RoleSnowflake
    ):
        state = StateService(ctx=ctx)
        snowflake_kwargs = {
            "channel_snowflake": int(ctx.channel.id),
            "guild_snowflake": int(ctx.guild.id),
            "member_snowflake": int(ctx.author.id),
        }
        do = DiscordObject(ctx=ctx)
        role_dict = await do.determine_from_target(target=role)
        pages = await AdministratorRoleService.toggle_administrator_role(
            role_dict=role_dict, snowflake_kwargs=snowflake_kwargs
        )
        await StateService.send_pages(title="Administrators", pages=pages, state=state)

    @commands.command(name="hero", help="Grant/revoke invincibility.")
    @guild_owner_predicator()
    @skip_help_discovery()
    async def invincibility_text_command(
        self,
        ctx: commands.Context,
        member: MemberSnowflake = commands.parameter(
            description="Tag a member or include their ID"
        ),
    ):
        state = StateService(ctx=ctx)
        do = DiscordObject(ctx=ctx)
        member_dict = await do.determine_from_target(target=member)
        where_kwargs = member_dict.get("columns", None)
        enabled = Invincibility.toggle_enabled()
        if enabled:
            Invincibility.add_invincible_member(**where_kwargs)
            await Invincibility.unrestrict(**where_kwargs)
            msg = (
                f"All moderation events have been forgiven "
                f"and invincibility has been enabled for {member_dict.get("mention", None)}."
            )
        else:
            Invincibility.remove_invincible_member(**where_kwargs)
            msg = f"Invincibility has been disabled for {member_dict.get("mention", None)}"
        return await state.end(success=msg)


async def setup(bot: DiscordBot):
    await bot.add_cog(GuildOwnerTextCommands(bot))
