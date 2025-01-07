''' hybrid.py The purpose of this program is to be an extension to a Discord bot to provide the command functionality from cd ..lucy./..lucy./.
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
from discord.utils import get
from discord import Embed
from discord.ext import commands, tasks
from PIL import Image
from random import randint
from typing import Optional
from ..utils.frames import extract_random_frames
from ..utils.add_watermark import add_watermark
from ..utils.average_score import average_score
from ..utils.combine import combine
from ..utils.create_completion import create_completion
from ..utils.draw_fingerprint import draw_fingerprint
from ..utils.draw_watermarked_molecule import draw_watermarked_molecule
from ..utils.get_mol import get_mol
from ..utils.google import google
from ..utils.gsrs import gsrs
from ..utils.script import script
from ..utils.unique_pairs import unique_pairs
from ..utils.tag import TagManager
#import aiomysql
import asyncio
import discord
#from googletrans import Translator, LANGUAGES
import io
import json
import os
import shlex
import traceback
from ..utils.helpers import *
import random

def at_home():
    async def predicate(ctx):
        return ctx.guild is not None and ctx.guild.id == 1300517536001036348
    return commands.check(predicate)

class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.indica = self.bot.get_cog('Indica')
        self.sativa = self.bot.get_cog('Sativa')
        self.tag_manager = TagManager(self.bot.db_pool)
        self.messages = []

    @commands.command(description='Change your role color using RGB values. Usage: between `!colorize 0 0 0` and `!colorize 255 255 255`')
    @at_home()
    async def colorize(self, ctx: commands.Context, r: Optional[str] = commands.parameter(default='149', description='Anything between 0 and 255.'), g: int = commands.parameter(default='165', description='Anything betwen 0 and 255.'), b: int = commands.parameter(default='165', description='Anything between 0 and 255.')):
        if not r.isnumeric():
            input_text_dict = {
                'type': 'text',
                'text': r
            }
            array = [
                {
                    'role': 'user',
                    'content': json.dumps(input_text_dict)
                }
            ]
            async for completion in create_completion(array):
                color_values = json.loads(completion)
                r = color_values['r']
                g = color_values['g']
                b = color_values['b']
        guildroles = await ctx.guild.fetch_roles()
        position = len(guildroles) - 12
        for arg in ctx.author.roles:
            if arg.name.isnumeric():
                await ctx.author.remove_roles(arg)
        for arg in guildroles:
            if arg.name.lower() == f'{r}{g}{b}':
                await ctx.author.add_roles(arg)
                await arg.edit(position=position)
                await ctx.send(f'I successfully changed your role color to {r}, {g}, {b}')
                return
        newrole = await ctx.guild.create_role(name=f'{r}{g}{b}', color=discord.Color.from_rgb(r, g, b), reason='new color')
        await newrole.edit(position=position)
        await ctx.author.add_roles(newrole)
        await ctx.send(f'I successfully changed your role color to {r}, {g}, {b}')

    @at_home()
    def get_language_code(self, language_name):
        language_name = language_name.lower()
        for lang_code, lang_name in LANGUAGES.items():
            if lang_name.lower() == language_name:
                return lang_code
        return None

#    @commands.command()
#    async def languages(self, ctx):
#        supported_languages = ', '.join(LANGUAGES.values())
#        await ctx.send(f'Supported languages are:\n{supported_languages}')
#
#    @commands.command()
#    async def translate(self, ctx, toggle: str, target_lang: str = 'english', source_lang: str = 'auto'):
#        if toggle.lower() == 'on':
#            target_lang_code = self.get_language_code(target_lang)
#            source_lang_code = self.get_language_code(source_lang)
#            if target_lang_code is None or source_lang_code is None:
#                await ctx.send(f'{ctx.author.mention}, please specify valid language names.')
#                return
#            self.user_translation_preferences[ctx.author.id] = (target_lang_code, source_lang_code)
#            await ctx.send(f'{ctx.author.mention}, translation enabled from {source_lang} to {target_lang}.')
#        elif toggle.lower() == 'off':
#            self.user_translation_preferences[ctx.author.id] = None
#            await ctx.send(f'{ctx.author.mention}, translation disabled.')
#        else:
#            await ctx.send(f'{ctx.author.mention}, please specify 'on' or 'off'.')

    @commands.hybrid_command(
        name='tag',
        description='Manage or retrieve tags. Sub-actions: add, update, remove, list, loop.'
    )
    async def tag_command(
        self,
        ctx: commands.Context,
        action: str = commands.parameter(
            default=None,
            description='Action: add, update, remove, list, or the name of a tag.'
        ),
        name: Optional[str] = commands.parameter(
            default=None,
            description='Name of the tag.'
        ),
        content: Optional[str] = commands.parameter(
            default=None,
            description='Content for the tag.'
        )
    ):
        attachment_url = ctx.message.attachments[0].url if ctx.message.attachments else None
        action = action.lower() if action else None

        if action == "add":
            if not name:
                await ctx.send("Usage: `!tag add <name> \"content\"`")
                return
            try:
                await self.tag_manager.add_tag(
                    name=name,
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id,
                    content=content,
                    attachment_url=attachment_url
                )
                await ctx.send(f"Tag `{name}` added successfully.")
            except Exception as e:
                logger.error(f"Error adding tag: {e}")
                await ctx.send("An error occurred while adding the tag.")

        elif action == "update":
            if not name or not content:
                await ctx.send("Usage: `!tag update <name> \"new content\"`")
                return
            try:
                result = await self.tag_manager.update_tag(
                    name=name,
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id,
                    content=content,
                    attachment_url=attachment_url
                )
                if result > 0:
                    await ctx.send(f"Tag `{name}` updated.")
                else:
                    await ctx.send(f"Tag `{name}` not found or you do not own it.")
            except Exception as e:
                logger.error(f"Error updating tag: {e}")
                await ctx.send("An error occurred while updating the tag.")

        elif action == "remove":
            if not name:
                await ctx.send("Usage: `!tag remove <name>`")
                return
            try:
                result = await self.tag_manager.delete_tag(
                    name=name,
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id
                )
                if result > 0:
                    await ctx.send(f"Tag `{name}` removed.")
                else:
                    await ctx.send(f"Tag `{name}` not found or you do not own it.")
            except Exception as e:
                logger.error(f"Error removing tag: {e}")
                await ctx.send("An error occurred while removing the tag.")

        elif action == "list":
            try:
                tags = await self.tag_manager.list_tags(ctx.guild.id, ctx.author.id)
                if not tags:
                    await ctx.send("You have no tags.")
                else:
                    tag_list = "\n".join(
                        f"**{t['name']}**: {t['content'] or t['attachment_url']}" for t in tags
                    )
                    await ctx.send(f"Your tags:\n{tag_list}")
            except Exception as e:
                logger.error(f"Error listing tags: {e}")
                await ctx.send("An error occurred while listing your tags.")

        else:
            try:
                tag = await self.tag_manager.get_tag(location_id=ctx.guild.id, name=action)
                tag_content = tag.get('content')
                tag_attachment_url = tag.get('attachment_url')

                if tag_content and tag_attachment_url:
                    async with aiohttp.ClientSession() as session:
                        async with session.get(tag_attachment_url) as response:
                            if response.status == 200:
                                file_bytes = await response.read()
                                filename = tag_attachment_url.split("/")[-1]
                                await ctx.send(
                                    content=tag_content,
                                    file=discord.File(io.BytesIO(file_bytes), filename=filename)
                                )
                            else:
                                await ctx.send(tag_content)
                elif tag_content:
                    await ctx.send(tag_content)
                elif tag_attachment_url:
                    async with aiohttp.ClientSession() as session:
                        async with session.get(tag_attachment_url) as response:
                            if response.status == 200:
                                file_bytes = await response.read()
                                filename = tag_attachment_url.split("/")[-1]
                                await ctx.send(file=discord.File(io.BytesIO(file_bytes), filename=filename))
                            else:
                                await ctx.send("Failed to retrieve the attachment.")
            except Exception as e:
                logger.error(f"Error retrieving tag: {e}")
                await ctx.send(f"Tag `{action}` not found.")

    @commands.command(name='script', description='Usage !script <NIV/ESV> <Book>.<Chapter>.<Verse>')
    @at_home()
    async def script(self, ctx: commands.Context, version: str, *, reference: str):
         try:
             await ctx.send(script(version, reference))
         except Exception as e:
             print(traceback.format_exc())

    @commands.command(name='load', hidden=True)
    @at_home()
    async def load(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')


    @commands.hybrid_command(name='draw', description='Usage: !draw glow <molecule> or !draw gsrs <molecule> or !draw shadow <molecule>.')
    @at_home()
    async def molecule(self, ctx: commands.Context, option: str = commands.parameter(default='glow', description='Compare `compare or Draw style `glow` `gsrs` `shadow`.'), *, molecules: str = commands.parameter(default=None, description='Any molecule'), quantity: int = commands.parameter(default=1, description='Quantity of glows')):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            if option == 'compare':
                if not molecules:
                    await ctx.send('No molecules provided.')
                    return
                args = shlex.split(molecules)
                pairs = unique_pairs(args)
                if not pairs:
                    embed = discord.Embed(description='No valid pairs found.')
                    await ctx.send(embed=embed)
                    return
                for pair in pairs:
                    mol = get_mol(pair[0])
                    refmol = get_mol(pair[1])
                    if mol is None or refmol is None:
                        embed = discord.Embed(description=f'One or both of the molecules {pair[0]} or {pair[1]} are invalid.')
                        await ctx.send(embed=embed)
                        continue
                    fingerprints = [
                        draw_fingerprint([mol, refmol]),
                        draw_fingerprint([refmol, mol])
                    ]
                    combined_image = combine(fingerprints, reversed(pair))
                    await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
            elif option == 'glow':
                if not molecules:
                    await ctx.send('No molecules provided.')
                    return
                args = shlex.split(molecules)
                fingerprints = []
                names = []
                molecule = get_mol(args[0])
                for _ in range(quantity.default):
                    names.append(args[0])
                    fingerprints.append(draw_fingerprint([molecule, molecule]))
                combined_image = combine(fingerprints, names)
                await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
            elif option == 'gsrs':
                if not molecules:
                    await ctx.send('No molecules provided.')
                    return
                args = shlex.split(molecules)
                for molecule_name in args:
                    if molecule_name is None:
                        await ctx.send(f'{molecule_name} is an unknown molecule.')
                        continue
                    watermarked_image = gsrs(molecule_name)
                    with io.BytesIO() as image_binary:
                        watermarked_image.save(image_binary, format='PNG')
                        image_binary.seek(0)
                        await ctx.send(file=discord.File(fp=image_binary, filename='watermarked_image.png'))
            elif option == 'shadow':
                if not molecules:
                    await ctx.send('No molecules provided.')
                    return
                args = shlex.split(molecules)
                mol = get_mol(args[0])
                if mol is None:
                    embed = discord.Embed(description='Invalid molecule name or structure.')
                    await ctx.send(embed=embed)
                    return
                image = draw_watermarked_molecule(mol)
                await ctx.send(file=discord.File(image, f'{args[0]}.png'))
            else:
                await ctx.send('Invalid option. Use `compare`, `glow`, `gsrs`, or `shadow`.')
        except Exception as e:
            logger.error(traceback.format_exc())
            await ctx.reply(e)

    @commands.hybrid_command(hidden=True)
    @at_home()
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.hybrid_command(name='search', description='Usage: !search <query>. Search Google.')
    @at_home()
    async def search(self, ctx: commands.Context, *, query: str = commands.parameter(default=None, description='Google search a query.')):
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
        results = google(query)
        embed = discord.Embed(title=f'Search Results for `{query}`', color=discord.Color.blue())
        for result in results:
            embed.add_field(name=result['title'], value=result['link'], inline=False)
        await ctx.send(embed=embed)

    @commands.hybrid_command(name='frame', description='')
    @at_home()
    async def frame(self, ctx: commands.Context):
        video_path = 'frogs.mov'
        output_dir = 'frames'
        frames = extract_random_frames(video_path, output_dir)
        for frame in frames:
            await ctx.send(file=discord.File(frame))

async def setup(bot: commands.Bot):
    await bot.add_cog(Hybrid(bot))
