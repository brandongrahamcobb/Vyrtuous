''' discord_message_service.py  The purpose of this program is to handle messages in Discord.

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
from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.inc.helpers import *

import asyncio
import discord

class DiscordMessageService:
    def __init__(self, bot: DiscordBot, db_pool):
        self.bot = bot
        
    async def send_message(self, ctx: commands.Context, *, content: str=None, file: discord.File=None, embed: discord.Embed=None, allowed_mentions: discord.AllowedMentions=discord.AllowedMentions.all()):
        can_send = ctx.guild and isinstance(ctx.channel, discord.abc.GuildChannel) and ctx.channel.permissions_for(ctx.guild.me).send_messages
        if can_send:
            try:
                await self._send_message(lambda **kw: ctx.reply(**kw), content=content, file=file, embed=embed, allowed_mentions=allowed_mentions)
            except discord.HTTPException as e:
                if e.code == 50035:
                    await self._send_message(lambda **kw: ctx.send(**kw), content=content, file=file, embed=embed, allowed_mentions=allowed_mentions)
                else:
                    raise
        else:
            await self.send_dm(ctx.author, content=content, file=file, embed=embed, allowed_mentions=allowed_mentions)
    
    async def _send_message(self, send_func, *, content: str=None, file: discord.File=None, embed: discord.Embed=None, allowed_mentions: discord.AllowedMentions=discord.AllowedMentions.all()):
        kwargs = {}
        if content:
            kwargs['content'] = content
        if file:
            kwargs['file'] = file
        if embed:
            kwargs['embed'] = embed
        if allowed_mentions is not None:
            kwargs['allowed_mentions'] = allowed_mentions
        await send_func(**kwargs)
         
class AppPaginator(discord.ui.View):
    def __init__(self, bot: DiscordBot, interaction: discord.Interaction, pages, *, timeout=60):
        super().__init__(timeout=timeout)
        self.bot = bot
        self.interaction = interaction
        self.pages = pages
        self.current_page = 0
        self.message = None
        
    async def start(self):
        if not self.pages:
            return await self.interaction.response.send_message(
                'There are no pages to display.',
                ephemeral=True
            )
        embed = self.pages[self.current_page]
        self.message = await self.interaction.response.send_message(
            embed=embed,
            ephemeral=True,
            view=self
        )
        
    @discord.ui.button(label='⬅️', style=discord.ButtonStyle.secondary)
    async def previous(self, interaction: discord.Interaction, button: discord.ui.Button):
        if interaction.user != self.interaction.user:
            return await interaction.response.send_message(
                "You cannot use this paginator.",
                ephemeral=True
            )
        if self.current_page > 0:
            self.current_page -= 1
            await interaction.response.edit_message(embed=self.pages[self.current_page], view=self)

    @discord.ui.button(label='➡️', style=discord.ButtonStyle.secondary)
    async def next(self, interaction: discord.Interaction, button: discord.ui.Button):
        if interaction.user != self.interaction.user:
            return await interaction.response.send_message(
                "You cannot use this paginator.",
                ephemeral=True
            )
        if self.current_page < len(self.pages) - 1:
            self.current_page += 1
            await interaction.response.edit_message(embed=self.pages[self.current_page], view=self)

    @discord.ui.button(label='⏹️', style=discord.ButtonStyle.danger)
    async def stop(self, interaction: discord.Interaction, button: discord.ui.Button):
        if interaction.user != self.interaction.user:
            return await interaction.response.send_message(
                "You cannot use this paginator.",
                ephemeral=True
            )
        await interaction.response.edit_message(view=None)
        self.stop()

class Paginator:
    def __init__(self, bot: DiscordBot, ctx_or_channel, pages):
        self.bot = bot
        self.ctx_or_channel = ctx_or_channel
        self.pages = pages
        self.current_page = 0
        self.message = None

    async def start(self):
        if isinstance(self.ctx_or_channel, commands.Context):
            messageable = self.ctx_or_channel
            author = self.ctx_or_channel.author
            author_only = True
        elif isinstance(self.ctx_or_channel, discord.abc.Messageable):
            messageable = self.ctx_or_channel
            author = None
            author_only = False
        else:
            raise TypeError(f"Expected Context or Messageable, got {type(self.ctx_or_channel)}")
        if not self.pages:
            await messageable.send('There are no pages to display.')
            return
        self.message = await messageable.send(embed=self.pages[self.current_page])
        await self.message.add_reaction('⬅️')
        await self.message.add_reaction('➡️')
        await self.message.add_reaction('⏹️')
        def check(reaction, user):
            if user.bot:
                return False
            if author_only:
                return (
                    user == author
                    and reaction.message.id == self.message.id
                    and str(reaction.emoji) in ['⬅️', '➡️', '⏹️']
                )
            return (
                reaction.message.id == self.message.id
                and str(reaction.emoji) in ['⬅️', '➡️', '⏹️']
            )
        while True:
            try:
                reaction, user = await self.bot.wait_for('reaction_add', timeout=60.0, check=check)
                if str(reaction.emoji) == '⬅️':
                    if self.current_page > 0:
                        self.current_page -= 1
                        await self.message.edit(embed=self.pages[self.current_page])
                elif str(reaction.emoji) == '➡️':
                    if self.current_page < len(self.pages) - 1:
                        self.current_page += 1
                        await self.message.edit(embed=self.pages[self.current_page])
                elif str(reaction.emoji) == '⏹️':
                    await self.message.clear_reactions()
                    break
                await self.message.remove_reaction(reaction.emoji, user)
            except asyncio.TimeoutError:
                await self.message.clear_reactions()
                break

