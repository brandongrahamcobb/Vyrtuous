"""!/bin/python3
list_overwrites.py The purpose of this program is to list channel overwrites.

Copyright (C) 2026  https://github.com/brandongrahamcobb/Vyrtuous.git

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

from vyrtuous.utils.messaging import emojis


def build_embed(obj) -> discord.Embed:
    member_count, role_count, total_count = 0, 0, 0
    obj_name = obj.name
    title = f"{emojis.get_random_emoji()} Overwrites for {obj_name}"
    if isinstance(obj, discord.Guild):
        for channel in obj.channels:
            for target_obj, overwrite in channel.overwrites.items():
                if any(v is not None for v in overwrite):
                    total_count += 1
                    if isinstance(target_obj, discord.Member):
                        member_count += 1
                    else:
                        role_count += 1
    elif isinstance(
        obj, (discord.TextChannel, discord.VoiceChannel, discord.StageChannel)
    ):
        for target_obj, overwrite in obj.overwrites.items():
            if any(v is not None for v in overwrite):
                total_count += 1
                if isinstance(target_obj, discord.Member):
                    member_count += 1
                else:
                    role_count += 1
    embed = discord.Embed(title=title, color=discord.Color.blue())
    embed.add_field(name="Role overwrites", value=str(role_count), inline=False)
    embed.add_field(name="Member overwrites", value=str(member_count), inline=False)
    embed.add_field(name="Total overwrites", value=str(total_count), inline=False)
    return embed
