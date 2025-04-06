''' hybrid.py The purpose of this program is to be an extension to a Discord bot to provide the command functionality to Vyrtuous.
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
from bs4 import BeautifulSoup
from collections import defaultdict
from discord.utils import get
from discord import Embed, File, app_commands
from discord.ext import commands, tasks
from lucy.utils.frames import extract_random_frames
from lucy.utils.add_watermark import add_watermark
from lucy.utils.average_score import average_score
from lucy.utils.combine import combine_gallery
from lucy.utils.create_batch_completion import BatchProcessor
from lucy.utils.create_completion import create_completion
from lucy.utils.draw_fingerprint import draw_fingerprint
from lucy.utils.draw_watermarked_molecule import draw_watermarked_molecule
from lucy.utils.game import Game
from lucy.utils.get_mol import construct_helm_from_peptide
from lucy.utils.get_mol import get_mol
from lucy.utils.get_mol import manual_helm_to_smiles
from lucy.utils.get_molecule_name import get_molecule_name
from lucy.utils.get_proximity import get_proximity
from lucy.utils.google import google
from lucy.utils.gsrs import gsrs
from lucy.utils.helpers import *
from lucy.utils.image import create_image, create_image_variation, edit_image
from lucy.utils.message import Message
from lucy.utils.paginator import Paginator
from lucy.utils.predicator import Predicator
from lucy.utils.script import script
from lucy.utils.stable_cascade import stable_cascade
from lucy.utils.tag import TagManager
from lucy.utils.unique_pairs import unique_pairs
from lucy.utils.usage import OpenAIUsageClient
from PIL import Image
from random import randint
from rdkit import Chem
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem import AllChem, Crippen
from random import choice
from typing import Dict, List, Optional
import asyncio
import datetime
import discord
from googletrans import Translator, LANGUAGES
import io
import json
import openai
import os
import pubchempy as pcp
import pytz
import re
import shlex
import time
import traceback
import uuid
from pyPept.sequence import Sequence, correct_pdb_atoms
from pyPept.molecule import Molecule
from pyPept.converter import Converter
from rdkit import Chem
import pubchempy as pcp
import re
import requests

class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = bot.conversations
        self.batch_processor = BatchProcessor(bot)
        self.game = Game(self.bot)
        self.predicator = Predicator(self.bot)
        self.tag_manager = TagManager(self.bot.db_pool)
        self.handler = Message(self.bot, self.config, self.conversations, self.bot.db_pool)
        self.loop_task: Optional[str] = None

    async def loop_tags(self, channel: discord.TextChannel):
        while True:
            try:
                loop_tags = await self.tag_manager.list_tags(channel.guild.id, tag_type='loop')
                if loop_tags:
                    random_tag = choice(loop_tags)
                    content_value = random_tag.get('content')
                    attachment_url_value = random_tag.get('attachment_url')
                    if content_value and attachment_url_value:
                        await channel.send(content_value)
                        await channel.send(attachment_url_value)
                    elif content_value:
                        await channel.send(content_value)
                    elif attachment_url_value:
                       await channel.send(attachment_url_value)
                await asyncio.sleep(300)
            except asyncio.CancelledError:
                break
            except Exception as e:
                logger.error(f'Error during loop_tags: {e}')

    def start_loop_task(self, channel: discord.TextChannel):
        if self.loop_task is None or self.loop_task.done():
            self.loop_task = asyncio.create_task(self.loop_tags(channel))

    def stop_loop_task(self):
        if self.loop_task and not self.loop_task.done():
            self.loop_task.cancel()
            self.loop_task = None


    def get_language_code(self, language_name):
        language_name = language_name.lower()
        for lang_code, lang_name in LANGUAGES.items():
            if lang_name.lower() == language_name:
                return lang_code
        return None

    @commands.hybrid_command(name="chat", description="Usage: chat <model> <prompt>")
    async def chat(
        self,
        ctx: commands.Context,
        model: str,
        prompt: str,
        new: bool = True,
        max_tokens: int = None,
        response_format: str = None,
        stop: str = None,
        store: bool = None,
        stream: bool = None,
        sys_input: str = None,
        temperature: float = None,
        top_p: float = None,
        use_history: bool = None,
        add_completion_to_history: bool = None,
    ):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            array = await self.handler.process_array(prompt, attachments=ctx.message.attachments)
            custom_id = f"{ctx.author.id}-{uuid.uuid4()}"
            request_data = {
                "completions": 1,
                "custom_id": custom_id,
                "input_array": array,
                "max_tokens": max_tokens if max_tokens is not None else OPENAI_MODEL_OUTPUT_LIMITS[model],
                "model": model,
                "response_format": response_format if response_format is not None else OPENAI_CHAT_RESPONSE_FORMAT,
                "stop": stop if stop is not None else self.config.get("openai_chat_stop", None),
                "store": store if store is not None else self.config.get("openai_chat_store", False),
                "stream": stream if stream is not None else self.config.get("openai_chat_stream", False),
                "sys_input": sys_input if sys_input is not None else self.config.get("openai_chat_sys_input", None),
                "temperature": temperature if temperature is not None else self.config.get("openai_chat_temperature", 0.7),
                "top_p": top_p if top_p is not None else self.config.get("openai_chat_top_p", 1.0),
                "use_history": use_history if use_history is not None else self.config.get("openai_chat_use_history", True),
                "add_completion_to_history": add_completion_to_history if add_completion_to_history is not None else self.config.get("openai_chat_add_completion_to_history", True),
            }
            if new:
                async for chat_completion in self.conversations.create_https_completion(**request_data):
                    if len(chat_completion) > 2000:
                        unique_filename = f'temp_{uuid.uuid4()}.txt'
                        with open(unique_filename, 'w') as f:
                            f.write(chat_completion)
                        await self.handler.send_message(ctx, content=None, file=discord.File(unique_filename))
                        os.remove(unique_filename)
                    else:
                        await self.handler.send_message(ctx, content=chat_completion)
            else:
                with open(PATH_OPENAI_REQUESTS, "a") as f:
                    f.write(json.dumps(request_data) + "\n")
                await self.handler.send_message(ctx, content="✅ Your request has been queued for weekend batch processing.")
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name='colorize', description=f'Usage: between `colorize 0 0 0` and `colorize 255 255 255` or `colorize <color>`')
    async def colorize(self, ctx: commands.Context, *, color: str = commands.parameter(default='blurple', description='Anything between 0 and 255 or a color.')):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        args = shlex.split(color)
        r = args[0]
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
        else:
            g = args[1]
            b = args[2]
        r = int(r)
        g = int(g)
        b = int(b)
        guildroles = await ctx.guild.fetch_roles()
        position = len(guildroles) - 12
        for arg in ctx.author.roles:
            if arg.name.isnumeric():
                await ctx.author.remove_roles(arg)
        for arg in guildroles:
            if arg.name.lower() == f'{r}{g}{b}':
                await ctx.author.add_roles(arg)
                await arg.edit(position=position)
                await self.handler.send_message(ctx, content=f'I successfully changed your role color to {r}, {g}, {b}')
                return
        newrole = await ctx.guild.create_role(name=f'{r}{g}{b}', color=discord.Color.from_rgb(r, g, b), reason='new color')
        await newrole.edit(position=position)
        await ctx.author.add_roles(newrole)
        await self.handler.send_message(ctx, content=f'I successfully changed your role color to {r}, {g}, {b}')

    @commands.hybrid_command(name='d', description=f'Usage: d 2 <mol> <mol> or d glow <mol> or d gsrs <mol> or d shadow <mole>.')
    async def d(
        self,
        ctx: commands.Context,
        option: str = commands.parameter(default='glow', description='Compare `compare or Draw style `glow` `gsrs` `shadow`.'),
        *,
        chems: str = commands.parameter(default=None, description='Any molecule'),
        dupes: int = commands.parameter(default=1, description='Quantity of glows'),
        backwards: bool = commands.parameter(default=False, description='Reverse'),
        linearity: bool = commands.parameter(default=False, description='Linearity'),
        rdkit_coords: bool = commands.parameter(default=True, description='rdDepictor'),
        r: int = commands.parameter(default=0, description='Rotation')
    ):
        if not self.predicator.is_release_mode_func(ctx):
            return
        rdkit_bool = rdkit_coords if isinstance(rdkit_coords, bool) else True
        quantity = dupes if isinstance(dupes, int) else 1
        rotation = r if isinstance(r, int) else 0
        reverse = backwards if isinstance(backwards, bool) else False
        if isinstance(quantity, commands.Parameter):
            quantity = 1
        if isinstance(rotation, commands.Parameter):
            rotation = 0
        async def function(linearity):
            async def get_attachment_text(ctx):
                if ctx.message.attachments:
                    content = await ctx.message.attachments[0].read()
                    return content.decode('utf-8')
                return None
            async def get_molecules(molecules):
                if not molecules:
                    await self.handler.send_message(ctx, content='No molecules provided.')
                    return
                name_match = re.search(r'"([^"]+)"$', molecules)
                if name_match:
                    name = name_match.group(1)
                    molecule_list = molecules[:name_match.start()].strip()
                else:
                    name = 'Untitled'
                    molecule_list = molecules
                if '.' in molecule_list:
                    molecule_parts = re.split(r'(?<!")\.(?!")', molecule_list)
                else:
                    molecule_parts = [molecule_list]
                fingerprints = []
                names = molecule_parts
                converted_smiles = []
                for mol in molecule_parts:
                    mol_obj = Chem.MolFromSmiles(mol)
                    if mol_obj:
                        smiles = mol
                    else:
                        compounds = pcp.get_compounds(mol, 'name')
                        if compounds:
                            smiles = compounds[0].isomeric_smiles
                        else:
                            helm = construct_helm_from_peptide(mol)
                            smiles = manual_helm_to_smiles(helm)
                    if not smiles:
                        embed = discord.Embed(description=f'Invalid molecule: {mol}')
                        await self.handler.send_message(ctx, content=None, file=None, embed=embed)
                        return
                    converted_smiles.append(smiles)
                molecule_objects = [get_mol(smiles, reverse=reverse) for smiles in converted_smiles]
                return molecule_objects, names, name
            if ctx.message.attachments:
                molecules = await get_attachment_text(ctx)
            molecules = chems or (await get_attachment_text(ctx) if ctx.message.attachments else None)
            if option == '2':
                molecule_objects, names, name = await get_molecules(molecules)
                if not molecules:
                    await self.handler.send_message(ctx, content='No molecules provided.')
                    return
                pairs = unique_pairs(names)
                if not pairs:
                    embed = discord.Embed(description='No valid pairs found.')
                    await self.handler.send_message(ctx, content=None, file=None, embed=embed)
                    return
                for pair in pairs:
                    mol = get_mol(pair[0], False)
                    refmol = get_mol(pair[1], False)
                    if mol is None or refmol is None:
                        embed = discord.Embed(description=f'One or both of the molecules {pair[0]} or {pair[1]} are invalid.')
                        await self.handler.send_message(ctx, content=None, file=None, embed=embed)
                        continue
                    else:
                        fingerprints = [
                            draw_fingerprint([mol, refmol], rdkit_bool),
                            draw_fingerprint([refmol, mol], rdkit_bool, rotation)
                        ]
                    if len(fingerprints) == 2:
                        linearity = True
                    combined_image = combine_gallery(fingerprints, names, name, 1, linearity)
                    await self.handler.send_message(ctx, content=None, file=discord.File(combined_image, f'molecule_comparison.png'))
            elif option == 'glow':
                molecule_objects, names, name = await get_molecules(molecules)
                fingerprints = [draw_fingerprint([mol_obj, mol_obj], rdkit_bool, rotation) for mol_obj in molecule_objects]
                combined_image = combine_gallery(fingerprints, names, name, quantity, linearity)
                await self.handler.send_message(ctx, content=None, file=discord.File(combined_image, 'molecule_comparison.png'))
            elif option == 'gsrs':
                if not molecules:
                    await self.handler.send_message(ctx, content='No molecules provided.')
                    return
                args = shlex.split(molecules)
                for molecule_name in args:
                    if molecule_name is None:
                        await self.handler.send_message(ctx, content=f'{molecule_name} is an unknown molecule.')
                        continue
                    watermarked_image = gsrs(molecule_name)
                    with io.BytesIO() as image_binary:
                        watermarked_image.save(image_binary, format='PNG')
                        image_binary.seek(0)
                        await self.handler.send_message(ctx, content=None, file=discord.File(fp=image_binary, filename='watermarked_image.png'))
            elif option == 'shadow':
                molecule_objects, names, name = await get_molecules(molecules)
                molecule_images = [draw_watermarked_molecule(mol_obj, rdkit_bool) for mol_obj in molecule_objects]
                if len(molecule_images) == 2:
                    linearity = True  # Set linearity to True for this case
                combined_image = combine_gallery(molecule_images, names, name, quantity, linearity)
                await self.handler.send_message(ctx, content=None, file=discord.File(combined_image, 'molecule_comparison.png'))
            else:
                await self.handler.send_message(ctx, content='Invalid option. Use `compare`, `glow`, `gsrs`, or `shadow`.')
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function(linearity)
                else:
                    await function(linearity)
            else:
                async with ctx.typing():
                    await function(linearity)

    @commands.hybrid_command()
    async def faction(self, ctx, action: str, *, faction_name: str = None):
        if not self.predicator.is_release_mode_func(ctx):
            return
        faction_input = str(faction_name) if faction_name is not None else None  # ✅ Moved outside
        async def function():
            user_id = ctx.author.id
            act = action.lower()
            user = await self.game.get_user(user_id)
            if not user:
                await self.handler.send_message(ctx, content="You have not interacted with the bot yet. Earn XP first!")
                return
            if act == "create":
                if not faction_input:
                    await self.handler.send_message(ctx, content="You must specify a faction name!")
                    return
                if user["faction_name"]:
                    await self.handler.send_message(ctx, content="You are already in a faction!")
                    return
                existing_faction = await self.game.get_faction(faction_input)
                if existing_faction:
                    await self.handler.send_message(ctx, content=f"Faction **{faction_input}** already exists!")
                    return
                success = await self.game.create_faction(faction_input, user_id)
                if success:
                    await self.handler.send_message(ctx, content=f"Faction **{faction_input}** has been created! 🎉")
                else:
                    await self.handler.send_message(ctx, content="There was an error creating the faction.")
                return
            elif act == "join":
                if not faction_input:
                    await self.handler.send_message(ctx, content="You must specify a faction to join!")
                    return
                if user["faction_name"]:
                    await self.handler.send_message(ctx, content="You are already in a faction!")
                    return
                response = await self.game.join_faction(faction_input, user_id)
                await self.handler.send_message(ctx, content=f"{ctx.author.mention}, {response}")
                return
            elif act == "leave":
                current_faction = user.get("faction_name")
                if not current_faction:
                    await self.handler.send_message(ctx, content="You are not in a faction!")
                    return
                await self.game.leave_faction(user_id, current_faction)
                await self.handler.send_message(ctx, content="You have left your faction and are now factionless.")
                return
            elif act == "switch":
                if not faction_input:
                    await self.handler.send_message(ctx, content="You must specify a faction to switch to!")
                    return
                if user["faction_name"] == faction_input:
                    await self.handler.send_message(ctx, content="You are already in that faction!")
                    return
                new_faction = await self.game.get_faction(faction_input)
                if not new_faction:
                    await self.handler.send_message(ctx, content="The faction you want to switch to does not exist!")
                    return
                if user["faction_name"]:
                    await self.game.leave_faction(user_id, user["faction_name"])
                response = await self.game.join_faction(faction_input, user_id)
                await self.handler.send_message(ctx, content=f"{ctx.author.mention}, you have switched to faction **{faction_input}**.")
                return
            elif act == "info":
                faction_to_use = faction_input or user.get("faction_name")
                if not faction_to_use:
                    await self.handler.send_message(ctx, content="You are not in a faction, and no faction name was provided!")
                    return
                faction_data = await self.game.get_faction(faction_to_use)
                if not faction_data:
                    await self.handler.send_message(ctx, content="Faction not found!")
                    return
                members = await self.game.get_faction_members(faction_to_use)
                members_count = len(members)
                await self.handler.send_message(ctx, content=
                    f"**Faction: {faction_to_use}**\n"
                    f"🔹 Level: {faction_data['level']}\n"
                    f"🔹 XP: {faction_data['xp']}\n"
                    f"🔹 Members: {members_count}"
                )
                return
            elif act == "leaderboard":
                factions = await self.game.get_faction_leaderboard()
                if not factions:
                    await self.handler.send_message(ctx, content="No factions found.")
                    return
                leaderboard = "**🏆 Faction Leaderboard:**\n"
                for i, faction in enumerate(factions[:10], start=1):
                    leaderboard += f"{i}. **{faction['name']}** - Level {faction['level']}, XP: {faction['xp']}\n"
                await self.handler.send_message(ctx, content=leaderboard)
                return
            else:
                await self.handler.send_message(ctx, content="Invalid action! Use `create`, `join`, `switch`, `info`, or `leaderboard`.")
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name='frame', description='Sends a frame from a number of animal cruelty footage sources.')
    async def frame(self, ctx: commands.Context):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            await ctx.interaction.response.defer(ephemeral=True)
            video_path = 'frogs.mov'
            output_dir = 'frames'
            frames = extract_random_frames(video_path, output_dir)
            for frame in frames:
                await self.handler.send_message(ctx, content=None, file=discord.File(frame))
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name="imagine")
    async def imagine(self, ctx, *, prompt: str):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            try:
                if ctx.message.attachments:
                    image_attachment = ctx.message.attachments[0]
                    image_bytes = await image_attachment.read()
                    image_file = discord.File(io.BytesIO(image_bytes), filename="uploaded_image.png")
                    message = await self.handler.send_message(ctx, content="Choose what to do with the image:", file=image_file)
                    await message.add_reaction("✅")
                    await message.add_reaction("❌")
                    await message.add_reaction("🖼️")
                    await message.add_reaction("🔲")
                    def check(reaction, user):
                        return user == ctx.author and str(reaction.emoji) in ["✅", "❌", "🖼️", "🔲"]
                    reaction, user = await self.bot.wait_for("reaction_add", check=check)
                    if str(reaction.emoji) == "✅":
                        await self.handler.send_message(ctx, content="Please upload a mask for editing, or confirm to use the full image as the mask.")
                        mask_msg = await self.bot.wait_for("message", check=lambda m: m.author == ctx.author)
                        if mask_msg.content.lower() == "confirm":
                            mask_file = image_attachment
                        else:
                            mask_file = mask_msg.attachments[0]
                        edited_image = await edit_image(image_file, mask_file, prompt)
                        if isinstance(edited_image, discord.File):
                            await self.handler.send_message(ctx, content="Here is your edited image with the mask:", file=edited_image)
                        else:
                            await self.handler.send_message(ctx, content=f"Error editing image: {edited_image}")
                    elif str(reaction.emoji) == "❌":
                        await self.handler.send_message(ctx, content="Edit canceled.")
                    elif str(reaction.emoji) == "🖼️":
                        variation = await create_image_variation(image_file, prompt)
                        if isinstance(variation, discord.File):
                            await self.handler.send_message(ctx, content="Here is your image variation:", file=variation)
                        else:
                            await self.handler.send_message(ctx, content=f"Error creating variation: {variation}")
                    elif str(reaction.emoji) == "🔲":
                        await self.handler.send_message(ctx, content="Please upload a mask image to use for editing.")
                        mask_msg = await self.bot.wait_for("message", check=lambda m: m.author == ctx.author)
                        mask_file = mask_msg.attachments[0]
                        edited_image = await edit_image(image_file, mask_file, prompt)
                        if isinstance(edited_image, discord.File):
                            await self.handler.send_message(ctx, content="Here is your edited image with the mask:", file=edited_image)
                        else:
                            await self.handler.send_message(ctx, content=f"Error editing image with mask: {edited_image}")
                else:
                    image_file = await create_image(prompt)
                    if isinstance(image_file, discord.File):
                        await self.handler.send_message(ctx, content="Here is your generated image:", file=image_file)
                    else:
                        await self.handler.send_message(ctx, content=f"Error generating image: {image_file}")
            except openai.OpenAIError as e:
                await self.handler.send_message(ctx, e.http_status)
                await self.handler.send_message(ctx, e.error)
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name='logp')
    async def logp(self, ctx: commands.Context, *, molecules: str):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            args = shlex.split(molecules)
            for arg in args:
                compounds = pcp.get_compounds(arg, 'name')
                compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
                log_p = Crippen.MolLogP(mol)
                await self.handler.send_message(ctx, content=f'Your octanol:water coefficient is: {log_p}')
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name='pic')
    async def pic(self, ctx: commands.Context, *, argument: str):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            file = stable_cascade(argument)
            if isinstance(file, discord.File):
                await self.handler.send_message(ctx, content=None, file=file)
            else:
                await self.handler.send_message(ctx, content=f"Error generating image: {file}")
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.command(name='script', description=f'Usage: lscript <NIV/ESV/Quran> <Book>.<Chapter>.<Verse>')
    async def script(self, ctx: commands.Context, version: str, *, reference: str):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            async with ctx.typing():
                await self.handler.send_message(ctx, content=script(version, reference))
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name='search', description=f'Usage: lsearch <query>. Search Google.')
    async def search(self, ctx: commands.Context, *, query: str = commands.parameter(default=None, description='Google search a query.')):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            results = google(query)
            embed = discord.Embed(title=f'Search Results for \"{query}\"', color=discord.Color.blue())
            for result in results:
                title, link = result.get("title", "No Title"), result.get("link", "No Link")
                embed.add_field(name=title, value=link, inline=False)
            await self.handler.send_message(ctx, content=None, file=None, embed=embed)
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name='sim')
    async def sim(self, ctx: commands.Context, *, molecules: str):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            args = shlex.split(molecules)
            similarity = get_proximity(get_mol(args[0]), get_mol(args[1]))
            await self.handler.send_message(ctx, content=similarity)
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name='smiles')
    async def smiles(self, ctx: commands.Context, *, molecules: str, reverse: bool = True):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            args = shlex.split(molecules)
            output = []
            for arg in args:
                compounds = pcp.get_compounds(arg, 'name')
                if compounds:
                    smiles_str = compounds[0].isomeric_smiles
                    output.append(f"{arg}: {smiles_str}")
                    continue
                else:
                    mol = Chem.MolFromSmiles(arg)
                    if mol:
                        smiles_str = Chem.MolToSmiles(mol, canonical=False)
                        output.append(f"{arg}: {smiles_str}")
                        continue
                try:
                    helm = construct_helm_from_peptide(arg)
                    smiles_str = manual_helm_to_smiles(helm)
                    output.append(f"{arg}: {smiles_str}")
                except Exception as e:
                    output.append(f"{arg}: Failed to generate molecule.")
                try:
                     mol = get_mol(arg, reverse=reverse)
                except Exception as e:
                    output.append(f"{arg}: Failed to generate molecule.")
            result = "\n".join(output)
            if len(result) > 2000:
                with open(f"smiles_{ctx.author.name}.txt", "w") as f:
                    f.write(result)
                await self.handler.send_message(ctx, content=None, file=discord.File(f"smiles_{ctx.author.name}.txt"))
            else:
                await self.handler.send_message(ctx, content=f"SMILES:\n```\n{result}\n```")
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.hybrid_command(name='tag', description='Manage or retrieve tags. Sub-actions: add, borrow, list, loop, rename, remove, update')
    async def tag_command(
        self,
        ctx: commands.Context,
        action: str = commands.parameter(default=None, description='Action: add, update, remove, list, loop.'),
        name: Optional[str] = commands.parameter(default=None, description='Name of the tag.'),
        content: Optional[str] = commands.parameter(default=None, description='Content for the tag (if applicable).'),
        tag_type: Optional[str] = commands.parameter(default=None, description='Optional tag type: default or loop.')
    ):
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            attachment_url = ctx.message.attachments[0].url if ctx.message.attachments else None
            act = action.lower() if action else None
            if act == 'add':
                resolved_tag_type = 'loop' if (tag_type and tag_type.lower() == 'loop') else 'default'
                if not name:
                    return await self.handler.send_message(ctx, content=f'Usage: \"{self.bot.command_prefix}tag add <name> \"content\" [loop]`')
                try:
                    await self.tag_manager.add_tag(
                        name=name,
                        location_id=ctx.guild.id,
                        owner_id=ctx.author.id,
                        content=content,
                        attachment_url=attachment_url,
                        tag_type=resolved_tag_type
                    )
                    await self.handler.send_message(ctx, content=f'Tag \"{name}\" (type: {resolved_tag_type}) added successfully.')
                except ValueError as ve:
                    await self.handler.send_message(ctx, content=str(ve))
                except Exception as e:
                    logger.error(f'Error adding tag: {e}')
                    await self.handler.send_message(ctx, content='An error occurred while adding the tag.')
            if act == "borrow":
                """
                Usage:
                    !tag borrow <tag_name>
                    Optionally, specify the original owner: !tag borrow <tag_name> @UserName
                """
                if not name:
                    return await self.handler.send_message(ctx, content=
                        f"Usage: `{self.bot.command_prefix}tag borrow <tag_name> [@owner]`"
                    )
                mentioned_users = ctx.message.mentions
                if mentioned_users:
                    owner = mentioned_users[0]
                    owner_id = owner.id
                else:
                    owner = None
                    owner_id = None
                try:
                    await self.tag_manager.borrow_tag(
                        tag_name=name,
                        location_id=ctx.guild.id,
                        borrower_id=ctx.author.id,
                        owner_id=owner_id,
                    )
                    if owner:
                        owner_display = owner.display_name
                        await self.handler.send_message(ctx, content=
                            f'You have successfully borrowed the tag "{name}" from {owner_display}.'
                        )
                    else:
                        await self.handler.send_message(ctx, content=
                            f'You have successfully borrowed the tag "{name}".'
                        )
                except ValueError as ve:
                    await self.handler.send_message(ctx, content=str(ve))
                except RuntimeError as re:
                    await self.handler.send_message(ctx, content=str(re))
                except Exception as e:
                    logger.error(f"Unexpected error during tag borrowing: {e}")
                    await self.handler.send_message(ctx, content=
                        "An unexpected error occurred while borrowing the tag."
                    )
            elif act == 'list':
                filter_tag_type = name.lower() if name and name.lower() in ('loop', 'default') else None
                try:
                    tags = await self.tag_manager.list_tags(
                        location_id=ctx.guild.id,
                        owner_id=ctx.author.id,
                        tag_type=filter_tag_type
                    )
                    if not tags:
                        await self.handler.send_message(ctx, content='No tags found.')
                    else:
                        tag_list = '\n'.join(f'**{t["name"]}**' for t in tags)
                        await self.handler.send_message(ctx, content=f'Tags:\n{tag_list}')
                except Exception as e:
                    logger.error(f'Error listing tags: {e}')
                    await self.handler.send_message(ctx, content='An error occurred while listing your tags.')
            elif act == 'remove':
                if not name:
                    return await self.handler.send_message(ctx, content=f'Usage: \"{self.bot.command_prefix}tag remove <name>`')
                try:
                    result = await self.tag_manager.delete_tag(
                        name=name,
                        location_id=ctx.guild.id,
                        owner_id=ctx.author.id
                    )
                    if result > 0:
                        await self.handler.send_message(ctx, content=f'Tag \"{name}\" removed.')
                    else:
                        await self.handler.send_message(ctx, content=f'Tag \"{name}\" not found or you do not own it.')
                except Exception as e:
                    logger.error(f'Error removing tag: {e}')
                    await self.handler.send_message(ctx, content='An error occurred while removing the tag.')
            elif act == 'loop':
                 if not name:
                     return await self.handler.send_message(ctx, content=f'Usage: \"{self.bot.command_prefix}tag loop on <#channel>` or \"{self.bot.command_prefix}tag loop off`')
                 if name.lower() == 'on':
                     channel = ctx.channel
                     if content and content.startswith('<#') and content.endswith('>'):
                         channel_id = int(content.strip('<#>'))
                         maybe_chan = self.bot.get_channel(channel_id)
                         if maybe_chan is not None:
                             channel = maybe_chan
                     try:
                         await self.tag_manager.set_loop_config(ctx.guild.id, channel.id, True)
                         self.start_loop_task(channel)
                         await self.handler.send_message(ctx, content=f'Looping enabled in {channel.mention}.')
                     except Exception as e:
                         logger.error(f'Error enabling loop: {e}')
                         await self.handler.send_message(ctx, content='Could not enable loop.')
                 elif name.lower() == 'off':
                     try:
                         await self.tag_manager.set_loop_config(ctx.guild.id, None, False)
                         self.stop_loop_task()
                         await self.handler.send_message(ctx, content='Looping disabled.')
                     except Exception as e:
                         logger.error(f'Error disabling loop: {e}')
                         await self.handler.send_message(ctx, content='Could not disable loop.')
            else:
                try:
                    tag = await self.tag_manager.get_tag(ctx.guild.id, act)
                    if tag:
                        content_value = tag.get('content')
                        attachment_url_value = tag.get('attachment_url')
                        if content_value and attachment_url_value:
                            await self.handler.send_message(ctx, content=content_value)
                            await self.handler.send_message(ctx, content=attachment_url_value)
                        elif content_value:
                            await self.handler.send_message(ctx, content=content_value)
                        elif attachment_url_value:
                            await self.handler.send_message(ctx, content=attachment_url_value)
                        else:
                            await self.handler.send_message(ctx, content=f'Tag \"{act}\" has no content.')
                    else:
                        await self.handler.send_message(ctx, content=f'Tag \"{act}\" not found.')
                except Exception as e:
                    logger.error(f'Error fetching tag \"{act}\": {e}')
                    await self.handler.send_message(ctx, content=f'An error occurred while fetching tag \"{act}\".')
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    @commands.command(name='wipe', description=f'Usage: lwipe <all|bot|commands|text|user>')
    @commands.has_permissions(manage_messages=True)
    async def wipe(self, ctx, option: str = None, limit: int = 100):
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        async def function():
            if limit <= 0 or limit > 100:
                return await self.handler.send_message(ctx, content='Limit must be between 1 and 100.')
            check_function = None
            if option == 'bot':
                check_function = lambda m: m.author == self.bot.user
            elif option == 'all':
                check_function = lambda m: True
            elif option == 'user':
                user = ctx.message.mentions[0] if ctx.message.mentions else None
                if user:
                    check_function = lambda m: m.author == user
                else:
                    return await self.handler.send_message(ctx, content='Please mention a user.')
            elif option == 'commands':
                check_function = lambda m: m.content.startswith(ctx.prefix)
            elif option == 'text':
                await self.handler.send_message(ctx, content='Provide text to delete messages containing it.')
                try:
                    msg_text = await self.bot.wait_for('message', timeout=30.0, check=lambda m: m.author == ctx.author)
                    check_function = lambda m: msg_text.content in m.content
                except asyncio.TimeoutError:
                    return await self.handler.send_message(ctx, content='You took too long to provide text. Cancelling operation.')
            else:
                return await self.handler.send_message(ctx, content='Invalid option.')
            total_deleted = 0
            while total_deleted < limit:
                deleted = await ctx.channel.purge(limit=min(limit - total_deleted, 10), check=check_function)
                if not deleted:
                    break
                total_deleted += len(deleted)
                await asyncio.sleep(1)
            if total_deleted > 0:
                await self.handler.send_message(ctx, content=f'Deleted {total_deleted} messages.')
            else:
                await self.handler.send_message(ctx, content='No messages matched the criteria.')
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
            await function()
        else:
            if ctx.channel and isinstance(ctx.channel, discord.abc.GuildChannel):
                permissions = ctx.channel.permissions_for(ctx.guild.me)
                if permissions.send_messages:
                    async with ctx.typing():
                        await function()
                else:
                    await function()
            else:
                async with ctx.typing():
                    await function()

    async def translate(self, ctx, toggle: str, target_lang: str = 'english', source_lang: str = 'auto'):
        if toggle.lower() == 'on':
            target_lang_code = self.get_language_code(target_lang)
            source_lang_code = self.get_language_code(source_lang)
            if target_lang_code is None or source_lang_code is None:
                await self.handler.send_message(ctx, content=f'{ctx.author.mention}, please specify valid language names.')
                return
            self.user_translation_preferences[ctx.author.id] = (target_lang_code, source_lang_code)
            await self.handler.send_message(ctx, content=f'{ctx.author.mention}, translation enabled from {source_lang} to {target_lang}.')
        elif toggle.lower() == 'off':
            self.user_translation_preferences[ctx.author.id] = None
            await self.handler.send_message(ctx, content=f'{ctx.author.mention}, translation disabled.')
        else:
            await self.handler.send_message(ctx, content=f'{ctx.author.mention}, please specify \'on\' or \'off\'.')

async def setup(bot: commands.Bot):
    await bot.add_cog(Hybrid(bot))
