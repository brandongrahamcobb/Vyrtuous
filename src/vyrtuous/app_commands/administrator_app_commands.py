"""!/bin/python3

admin_app_commands.py A discord.py cog containing administrative commands for the Vyrtuous bot.

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

import discord
from discord import app_commands
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import at_home
from vyrtuous.listing import list_overwrites, list_server_mutes
from vyrtuous.models import multi_converter
from vyrtuous.utils.messaging import emojis
from vyrtuous.utils.messaging.tick import Tick
from vyrtuous.utils.moderation import (
    clear_service,
    server_mute_service,
    voice_mute_service,
)
from vyrtuous.utils.users import coordinator_service, moderator_service
from vyrtuous.view.cancel_confirm_view import VerifyView


class AdministratorAppCommands(commands.Cog):
    PERMISSION_LEVEL = "Administrator"

    def __init__(self, bot: DiscordBot):
        self.__bot = bot

    async def interaction_check(self, interaction: discord.Interaction):
        if interaction.guild is None:
            raise commands.CheckFailure(
                "This command must be executed inside a server."
            )
        if interaction.channel is None:
            raise commands.CheckFailure(
                "This command must be executed in a valid channel."
            )
        await moderator_service.check_minimum_role(
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=interaction.user.id,
            lowest_role=self.PERMISSION_LEVEL,
        )
        return True

    @app_commands.command(name="clear", description="Reset records.")
    @app_commands.describe(
        category="Specify one of: `admin`, `alias`, `all`, `automute`, `ban`, `coord`, `flag`, `mod`, `tmute`, `stream` or `vmute`.",
        scope="Specify one of: `auto`, `server`, `user`",
        target="Specify one of: channel ID/mention, member ID/mention, or server ID.",
    )
    async def clear_channel_access_app_command(
        self,
        interaction: discord.Interaction,
        target: str,
        category: str | None = "all",
        scope: str | None = "user",
    ):
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        if target == "all":
            obj = None
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        view = VerifyView(
            author_snowflake=interaction.user.id,
            category=str(category),
            obj=obj,
        )
        embed = view.build_embed()
        await tick.end(success=embed, view=view)
        await view.wait()
        tick = Tick(bot=self.__bot, interaction=interaction)
        message = await interaction.original_response()
        message_snowflake = message.id
        msg = await clear_service.clear(
            author_snowflake=interaction.user.id,
            category=str(category),
            guild_snowflake=interaction.guild.id,
            message_snowflake=message_snowflake,
            message_channel_snowflake=message.channel.id,
            obj=obj,
            target=scope or "user",
            view=view,
        )
        return await tick.end(success=msg)

    # TODO: Make guild agnostic
    @app_commands.command(name="coord", description="Grant/revoke coords.")
    @app_commands.describe(
        channel="Specify a channel ID/mention",
        member="Specify a member ID/mention",
    )
    async def toggle_coordinator_app_command(
        self,
        interaction: discord.Interaction,
        member: str,
        channel: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        target_member = multi_converter.transform(
            interaction=interaction, argument=member
        )
        if isinstance(target_member, int):
            member_snowflake = target_member
        elif isinstance(target_member, discord.Member):
            member_snowflake = target_member.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        target_channel = multi_converter.transform(
            interaction=interaction, argument=channel
        )
        if target_channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = interaction.channel.id
        elif isinstance(
            target_channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = target_channel.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        await moderator_service.has_equal_or_lower_role(
            target_member_snowflake=member_snowflake,
            member_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=interaction.guild.id,
        )
        message = await interaction.original_response()
        message_snowflake = message.id
        msg = await coordinator_service.toggle_coordinator(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=interaction.guild.id,
            member_snowflake=member_snowflake,
            message_snowflake=message_snowflake,
            message_channel_snowflake=channel_snowflake,
        )
        return await tick.end(success=msg)

    @app_commands.command(name="ow", description="Overwrite stats.")
    @app_commands.describe(target="Specify one of: a channel ID/mention or server ID.")
    async def list_overwrites_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        obj = multi_converter.transform(interaction=interaction, argument=target)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.channel
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        embed = list_overwrites.build_embed(obj=obj)
        return await tick.end(success=embed)

    @app_commands.command(name="rmute", description="Room mute (except yourself).")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
        duration="Specify a duration in m/h/d.",
        reason="Specify a reason.",
    )
    async def channel_mute_app_command(
        self,
        interaction: discord.Interaction,
        channel: str | None,
        duration: str | None = "1h",
        reason: str | None = "No reason provided.",
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        target_channel = multi_converter.transform(
            interaction=interaction, argument=channel
        )
        if target_channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = interaction.channel.id
        elif isinstance(
            target_channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = target_channel.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        pages = await voice_mute_service.channel_mute(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            duration_value=duration or "1h",
            guild_snowflake=interaction.guild.id,
            reason=reason or "No reason provided.",
        )
        return await tick.end(success=pages)

    @app_commands.command(name="rmv", description="VC move.")
    @app_commands.describe(
        source_channel="Specify a `from` channel ID/mention.",
        target_channel="Specify a `to` channel ID/mention.",
    )
    async def channel_move_all_app_command(
        self,
        interaction: discord.Interaction,
        source_channel: str,
        target_channel: str,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(
                warning="This command must be executed within a server."
            )
        target = multi_converter.transform(
            interaction=interaction, argument=target_channel
        )
        if target is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            target_channel_obj = interaction.channel
        elif isinstance(
            target,
            (discord.VoiceChannel, discord.StageChannel),
        ):
            target_channel_obj = target
        else:
            return await tick.end(warning="This command must target a valid channel.")
        source = multi_converter.transform(
            interaction=interaction, argument=source_channel
        )
        if source is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            source_channel_obj = interaction.channel
        elif isinstance(
            source,
            (discord.VoiceChannel, discord.StageChannel),
        ):
            source_channel_obj = source
        else:
            return await tick.end(warning="This command must target a valid channel.")
        if isinstance(
            target_channel_obj, (discord.VoiceChannel, discord.StageChannel)
        ) and isinstance(
            source_channel_obj, (discord.VoiceChannel, discord.StageChannel)
        ):
            failed, moved = [], []
            for member in source_channel_obj.members:
                try:
                    await member.move_to(target_channel_obj)
                    moved.append(member)
                except discord.Forbidden as e:
                    failed.append(member)
                    self.__bot.logger.warning(
                        f"Unable to move member "
                        f"{member.display_name} ({member.id}) from channel "
                        f"{source_channel_obj.name} ({source_channel_obj.id}) to channel "
                        f"{target_channel_obj.name} ({target_channel_obj.id}) in guild "
                        f"{interaction.guild.name} ({interaction.guild.id}). "
                        f"{str(e).capitalize()}"
                    )
            embed = discord.Embed(
                title=f"{emojis.get_random_emoji()} "
                f"Moved {source_channel_obj.mention} to "
                f"{target_channel_obj.mention}",
                color=discord.Color.green(),
            )
            if moved:
                embed.add_field(
                    name=f"Successfully Moved (`{len(moved)}`)",
                    value="\n".join(member.mention for member in moved),
                    inline=False,
                )
            else:
                embed.add_field(name="Successfully Moved", value="None", inline=False)
            if failed:
                embed.add_field(
                    name=f"Failed to Move ({len(failed)})",
                    value="\n".join(member.mention for member in failed),
                    inline=False,
                )
            embed.set_footer(
                text=f"Moved from {source_channel_obj.name} to {target_channel_obj.name}"
            )
            return await tick.end(success=embed)
        return await tick.end(warning="No members moved.")

    @app_commands.command(name="smute", description="Server mute/server unmute.")
    @app_commands.describe(
        member="Specify a member ID/mention.", reason="Specify a reason."
    )
    async def toggle_server_mute_app_command(
        self,
        interaction: discord.Interaction,
        member: str,
        reason: str | None = "No reason provided.",
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a server.")
        target = multi_converter.transform(interaction=interaction, argument=member)
        if isinstance(target, int):
            member_snowflake = target
        elif isinstance(target, discord.Member):
            member_snowflake = target.id
        else:
            return await tick.end(warning=f"Member {member} was not found.")
        if not isinstance(
            interaction.channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            return await tick.end(
                warning="This command must be executed in a valid channel."
            )
        msg = await server_mute_service.toggle_server_mute(
            author_snowflake=interaction.user.id,
            channel_snowflake=interaction.channel.id,
            guild_snowflake=interaction.guild.id,
            member_snowflake=member_snowflake,
            reason=reason or "No reason provided.",
        )
        return await tick.end(success=msg)

    @app_commands.command(name="smutes", description="List mutes.")
    @app_commands.describe(target="Specify one of: a member ID/mention or a server ID.")
    async def list_server_mutes_app_command(
        self,
        interaction: discord.Interaction,
        target: str | None,
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if target == "all":
            obj = None
        elif target is None:
            obj = interaction.guild
        else:
            obj = multi_converter.transform(interaction=interaction, argument=target)
        is_at_home = at_home(source=interaction)
        pages = await list_server_mutes.build_pages(obj=obj, is_at_home=is_at_home)
        return await tick.end(success=pages)

    @app_commands.command(name="xrmute", description="Unmute all.")
    @app_commands.describe(
        channel="Specify a channel ID/mention.",
    )
    async def channel_unmute_app_command(
        self, interaction: discord.Interaction, channel: str | None
    ) -> discord.Message:
        tick = Tick(bot=self.__bot, interaction=interaction)
        if interaction.guild is None:
            return await tick.end(warning="This command must be executed in a guild.")
        target_channel = multi_converter.transform(
            interaction=interaction, argument=channel
        )
        if target_channel is None:
            if interaction.channel is None:
                return await tick.end(
                    warning="This command must target a valid channel."
                )
            channel_snowflake = interaction.channel.id
        elif isinstance(
            target_channel,
            (discord.VoiceChannel, discord.StageChannel, discord.TextChannel),
        ):
            channel_snowflake = target_channel.id
        else:
            return await tick.end(warning="This command must target a valid channel.")
        pages = await voice_mute_service.channel_unmute(
            author_snowflake=interaction.user.id,
            channel_snowflake=channel_snowflake,
            guild_snowflake=interaction.guild.id,
            target="user",
        )
        return await tick.end(success=pages)


async def setup(bot: DiscordBot):
    await bot.add_cog(AdministratorAppCommands(bot=bot))
