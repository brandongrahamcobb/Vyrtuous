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
from discord import Embed, app_commands
from discord.ext import commands, tasks
from PIL import Image
from random import randint
from typing import Optional
from lucy.utils.frames import extract_random_frames
from lucy.utils.add_watermark import add_watermark
from lucy.utils.average_score import average_score
from lucy.utils.combine import combine
from lucy.utils.create_batch_completion import BatchProcessor
from lucy.utils.create_completion import create_completion
from lucy.utils.draw_fingerprint import draw_fingerprint
from lucy.utils.draw_watermarked_molecule import draw_watermarked_molecule
from lucy.utils.game import Game
from lucy.utils.get_mol import get_mol
from lucy.utils.get_molecule_name import get_molecule_name
from lucy.utils.get_proximity import get_proximity
from lucy.utils.google import google
from lucy.utils.gsrs import gsrs
from lucy.utils.helpers import *
from lucy.utils.paginator import Paginator
from lucy.utils.predicator import Predicator
from lucy.utils.script import script
from lucy.utils.tag import TagManager
from lucy.utils.unique_pairs import unique_pairs
from lucy.utils.usage import OpenAIUsageClient
from rdkit import Chem
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem import AllChem
from random import choice
from typing import Dict, List, Optional
#import aiomysql
import asyncio
import datetime
import discord
#from googletrans import Translator, LANGUAGES
import io
import json
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
from lucy.utils.get_mol import construct_helm_from_peptide
from lucy.utils.get_mol import manual_helm_to_smiles
import pubchempy as pcp
import re
import requests

async def get_attachment_text(ctx):
    if ctx.message.attachments:
        content = await ctx.message.attachments[0].read()
        return content.decode('utf-8')
    return None

class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.conversations = bot.conversations
        self.batch_processor = BatchProcessor(bot)
        self.game = Game(self.bot)
        self.predicator = Predicator(self.bot)
        self.tag_manager = TagManager(self.bot.db_pool)
        self.messages = []
        self.stacks = {}
        self.creatine = get_mol('Creatine')
        self.theanine = get_mol('L-theanine')
        self.trimethylglycine = get_mol('Trimethylglycine')
        self.serotonin = get_mol('Serotonin')

#    def get_language_code(self, language_name):
#        language_name = language_name.lower()
#        for lang_code, lang_name in LANGUAGES.items():
#            if lang_name.lower() == language_name:
#                return lang_code
#        return None

    def get_user_stack(self, user_id: int):
        return self.stacks.get(user_id, [])

    def set_user_stack(self, user_id: int, molecule_names: List[str]):
        mol_objects = [get_mol(name) for name in molecule_names]
        self.stacks[user_id] = mol_objects

    @commands.hybrid_command(name="batch_results", with_app_command=True)
    async def batch_results(self, ctx: commands.Context):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        responses = self.batch_processor.get_user_responses(ctx.author)
        if responses:
            response_text = "\n\n".join(responses)
            if len(response_text) > 2000:
                await ctx.send("Responses are too long. Sending as a file.")
                with open(f"batch_{ctx.author.name}.txt", "w") as f:
                    f.write(response_text)
                await ctx.send(file=discord.File(f"batch_{ctx.author.name}.txt"))
            else:
                await ctx.send(response_text)
        else:
            await ctx.send("No batch responses available.")

    @commands.hybrid_command(name="chat", with_app_command=True)
    async def chat(
        self,
        ctx: commands.Context,
        model: str,
        prompt: str,
        new: bool = True,  # True = Execute now, False = Queue for batch processing
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
        """
        Slash command to submit an OpenAI chat request.
        - Users **must specify `model` and `prompt`**.
        - Other parameters are **optional** and default to values from `config` unless overridden.
        - **If `new=True`**, request executes immediately.
        - **If `new=False`**, request is queued for batch processing.
        """
        await ctx.defer()
        input_array = [{"role": "user", "content": prompt}]
        custom_id = f"{ctx.author.id}-{uuid.uuid4()}"  # Unique ID per request
        request_data = {
            "completions": 1,
            "custom_id": custom_id,
            "input_array": input_array,
            "max_tokens": max_tokens if max_tokens is not None else self.config.get("openai_chat_max_tokens", 512),
            "model": model,  # User-specified model
            "response_format": response_format if response_format is not None else self.config.get("openai_chat_response_format", "json"),
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
                    await ctx.send(file=discord.File(unique_filename))
                    os.remove(unique_filename)
                else:
                    await ctx.send(chat_completion)
        else:
            with open(PATH_OPENAI_REQUESTS, "a") as f:
                f.write(json.dumps(request_data) + "\n")
            await ctx.send("✅ Your request has been queued for weekend batch processing.")

    @commands.hybrid_command(name='colorize', description=f'Usage: between `lcolorize 0 0 0` and `lcolorize 255 255 255` or `l colorize <color>`')
    @commands.has_permissions(manage_roles=True) # Do you have manage_roles permissions?
    async def colorize(self, ctx: commands.Context, r: str = commands.parameter(default='blurple', description='Anything between 0 and 255 or a color.'), *, g: str = commands.parameter(default='147', description='Anything betwen 0 and 255.'), b: str = commands.parameter(default='165', description='Anything between 0 and 255.')):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
#        if not await self.predicator.is_at_home_func(ctx.guild.id):
 #           return
        if not self.predicator.is_release_mode_func(ctx):
            return
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
                await ctx.send(f'I successfully changed your role color to {r}, {g}, {b}')
                return
        newrole = await ctx.guild.create_role(name=f'{r}{g}{b}', color=discord.Color.from_rgb(r, g, b), reason='new color')
        await newrole.edit(position=position)
        await ctx.author.add_roles(newrole)
        await ctx.send(f'I successfully changed your role color to {r}, {g}, {b}')

    @commands.hybrid_command(name='d', description=f'Usage: ld 2 <mol> <mol> or ld glow <mol> or ld gsrs <mol> or ld shadow <mole>.')
    async def d(
        self,
        ctx: commands.Context,
        option: str = commands.parameter(default='glow', description='Compare `compare or Draw style `glow` `gsrs` `shadow`.'),
        *,
        molecules: str = commands.parameter(default=None, description='Any molecule'),
        quantity: int = commands.parameter(default=1, description='Quantity of glows'),
        reverse: bool = commands.parameter(default=False, description='Reverse')
    ):
        try:
            async def get_attachment_text(ctx):
                if ctx.message.attachments:
                    content = await ctx.message.attachments[0].read()
                    return content.decode('utf-8')
                return None
            if ctx.interaction:
                async with ctx.typing():
                    await ctx.interaction.response.defer(ephemeral=True)
#            if not await self.predicator.is_at_home_func(ctx.guild.id):
 #               return
            if not self.predicator.is_release_mode_func(ctx):
                return
            if ctx.message.attachments:
                molecules = await get_attachment_text(ctx)
            if option == '2':
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
                    mol = get_mol(pair[0], reverse)
                    refmol = get_mol(pair[1], reverse)
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
            
                # Split molecules by periods outside of quotes, but preserve periods inside quotes
                molecule_parts = re.split(r'(?<!")\.(?!")', molecules)  # Split by period but avoid splitting inside quotes
            
                # Extract the quoted name from the end of the input, including the period inside quotes
                name_match = re.search(r'"([^"]+)"$', molecules)
                
                if name_match:
                    name = name_match.group(1)  # Extract the name from the quotes
                    molecule_parts = molecule_parts[:-1]  # Remove the last part which is the name
                else:
                    name = 'Untitled'  # Default name if no quoted name is found
                
                # Initialize lists for storing results
                fingerprints = []
                names = [name]  # Add the extracted name to the list
            
                # Process the molecule list
                converted_smiles = []  # Store SMILES representations
                for mol in molecule_parts:
                    compounds = pcp.get_compounds(mol, 'name')  # Try to fetch from PubChem
                    mol_obj = Chem.MolFromSmiles(mol)
                    if mol_obj:
                        smiles = mol  # Valid SMILES provided directly
                    else:
                        if compounds:
                            smiles = compounds[0].isomeric_smiles  # Get the SMILES representation
                        else:
                            helm = construct_helm_from_peptide(mol)
                            smiles = manual_helm_to_smiles(helm)  # Fallback for peptides
            
                    if not smiles:
                        embed = discord.Embed(description=f'Invalid molecule: {mol}')
                        await ctx.send(embed=embed)
                        return
                    
                    converted_smiles.append(smiles)
            
                full_smiles = ".".join(converted_smiles)
                smiles_comparison = [full_smiles, full_smiles]
                molecule_objects = [get_mol(full_smiles, reverse=reverse) for _ in range(2)]
                fingerprints.append(draw_fingerprint(molecule_objects))
               
                # Combine the image and send it
                combined_image = combine(fingerprints, names)
                await ctx.send(file=discord.File(combined_image, 'molecule_comparison.png'))
#                if not molecules:
#                    await ctx.send('No molecules provided.')
#                    return
#                args = shlex.split(molecules)[:-1]
#                fingerprints = []
#                names = []
#                smiles_list = molecules.split('.')  # Split input by "."
#                converted_smiles = []  # Store SMILES representations
#                for mol in smiles_list:
#                    compounds = pcp.get_compounds(mol, 'name')  # Try to fetch from PubChem
#                    mol_obj = Chem.MolFromSmiles(mol)
#                    if mol_obj:
#                        smiles = mol  # Valid SMILES provided directly
#                    else:
#                        if compounds:
#                            smiles = compounds[0].isomeric_smiles  # Get the SMILES representation
#                        else:
#                            helm = construct_helm_from_peptide(mol)
#                            smiles = manual_helm_to_smiles(helm)  # Fallback for peptides
#                    if not smiles:
#                        embed = discord.Embed(description=f'Invalid molecule: {mol}')
#                        await ctx.send(embed=embed)
#                        return
#                    converted_smiles.append(smiles)
#                if len(smiles_list) == 1:
#                    names.append(mol)
#                else:
#                    names.append(molecules[-1])
##                    names.append(mol)
#                full_smiles = ".".join(converted_smiles)
#                smiles_comparison = [full_smiles, full_smiles]
#                molecule_objects = [get_mol(full_smiles, reverse=reverse) for _ in range(2)]
#                fingerprints.append(draw_fingerprint(molecule_objects))
#                combined_image = combine(fingerprints, names)
#                await ctx.send(file=discord.File(combined_image, 'molecule_comparison.png'))
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
                mol = get_mol(args[0], reverse)
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

    @commands.hybrid_command(name='frame', description='Sends a frame from a number of animal cruelty footage sources.')
    async def frame(self, ctx: commands.Context):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
#        if not await self.predicator.is_at_home_func(ctx.guild.id):
 #           return
        if not self.predicator.is_release_mode_func(ctx):
            return
        video_path = 'frogs.mov'
        output_dir = 'frames'
        frames = extract_random_frames(video_path, output_dir)
        for frame in frames:
            await ctx.send(file=discord.File(frame))

    @commands.command()
    async def level(self, ctx, member: discord.Member = None):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        self.game.load_users()
        member = member or ctx.author
        user_id = int(member.id)

        if user_id in self.game.users:
            level = self.game.users[user_id]["level"]
            xp = self.game.users[user_id]["xp"]
            xp_needed_for_next_level = self.game.get_xp_for_level(level + 1) - xp
            await ctx.send(f"{member.mention} is at level {level} with {xp:.2f} XP.\n"
                           f"You need {xp_needed_for_next_level:.2f} XP to reach level {level + 1}.")
        else:
            await ctx.send(f"{member.mention} has not interacted with the bot yet.")

    @commands.command()
    async def leaderboard(self, ctx):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        self.game.load_users()
        # Sort users by level, then by XP
        sorted_users = sorted(self.game.users.items(), key=lambda item: (item[1]["level"], item[1]["xp"]), reverse=True)
        
        # Create leaderboard message
        leaderboard_message = "Leaderboard:\n"
        for i, (user_id, data) in enumerate(sorted_users[:10], start=1):
            member = ctx.guild.get_member(int(user_id))
            member_name = member.display_name if member else "Unknown User"
            level = data["level"]
            xp = data["xp"]
            leaderboard_message += f"{i}. {member_name} - Level {level}, {xp:.2f} XP\n"

        await ctx.send(leaderboard_message)

    @commands.hybrid_command(name='logp')
    async def logp(self, ctx: commands.Context, *, molecules: str):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        args = shlex.split(molecules)
        for arg in args:
            compounds = pcp.get_compounds(arg, 'name')
            compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
            mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
            log_p = Crippen.MolLogP(mol)
            await ctx.send(f'Your octanol:water coefficient is: {log_p}')

    @commands.command(name='script', description=f'Usage: lscript <NIV/ESV/Quran> <Book>.<Chapter>.<Verse>')
    async def script(self, ctx: commands.Context, version: str, *, reference: str):
        try:
            if ctx.interaction:
                async with ctx.typing():
                    await ctx.interaction.response.defer(ephemeral=True)
#            if not await self.predicator.is_at_home_func(ctx.guild.id):
 #               return
            if not self.predicator.is_release_mode_func(ctx):
                return
            await ctx.send(script(version, reference))
        except Exception as e:
            print(traceback.format_exc())

    @commands.hybrid_command(name='search', description=f'Usage: lsearch <query>. Search Google.')
    async def search(self, ctx: commands.Context, *, query: str = commands.parameter(default=None, description='Google search a query.')):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
#        if not await self.predicator.is_at_home_func(ctx.guild.id):
 #           return
  #      if not self.predicator.is_release_mode_func(ctx):
   #         return
        results = google(query)
        embed = discord.Embed(title=f'Search Results for \"{query}\"', color=discord.Color.blue())
        for result in results:
            title, link = result.get("title", "No Title"), result.get("link", "No Link")
            embed.add_field(name=title, value=link, inline=False)
        await ctx.send(embed=embed)

    @commands.hybrid_command(name='sim')
    async def sim(self, ctx: commands.Context, *, molecules: str):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        args = shlex.split(molecules)
        similarity = get_proximity(get_mol(args[0]), get_mol(args[1]))
        await ctx.send(similarity)

    @commands.hybrid_command(name='smiles')
    async def smiles(self, ctx: commands.Context, *, molecules: str, reverse: bool = True):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        try:
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
                await ctx.send(file=discord.File(f"smiles_{ctx.author.name}.txt"))
            else:
                await ctx.send(f"SMILES:\n```\n{result}\n```")
        except Exception as e:
            await ctx.send(f"Error: {e}")
#                compounds = pcp.get_compounds(arg, 'name')
#                if compounds:
#                    compound = compounds[0]
#                    smiles_str = compound.isomeric_smiles
#                    output.append(f"{arg}: {smiles_str}")
#                    continue
#                else:
#                    mol = Chem.MolFromSmiles(arg)
#                    if mol:
#                        smiles_str = Chem.MolToSmiles(mol)
#                        output.append(f"{arg}: {smiles_str}")
#                        continue
#                helm = construct_helm(arg, reverse=reverse)
#                try:
#                    smiles_str = helm_to_smiles_manual(helm)
#                    output.append(f"{arg}: {smiles_str}")
#                except Exception as e:
#                    logger.warning(f"Direct molecule conversion failed for '{arg}': {e}")
#        except Exception as e:
#            logger.warning(f"Direct molecule conversion failed for '{arg}': {e}")
                # Otherwise, assume it's a peptide.
#        mol = get_mol(arg)
#        if mol is None:
#            output.append(f"{arg}: Failed to generate molecule.")
#        else:
#            smiles_str = Chem.MolToSmiles(mol)
#            output.append(f"{arg}: {smiles_str}")
#        result = "\n".join(output)
#        else:
#            await ctx.send(f"SMILES:\n```\n{result}\n```")
##    @commands.command()
#    async def languages(self, ctx):
#        supported_languages = ', '.join(LANGUAGES.values())
#        await ctx.send(f'Supported languages are:\n{supported_languages}')
#

    @commands.hybrid_command(name='tag', description='Manage or retrieve tags. Sub-actions: add, borrow, list, loop, rename, remove, update')
    async def tag_command(
        self,
        ctx: commands.Context,
        action: str = commands.parameter(default=None, description='Action: add, update, remove, list, loop.'),
        name: Optional[str] = commands.parameter(default=None, description='Name of the tag.'),
        content: Optional[str] = commands.parameter(default=None, description='Content for the tag (if applicable).'),
        tag_type: Optional[str] = commands.parameter(default=None, description='Optional tag type: default or loop.')
    ):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        attachment_url = ctx.message.attachments[0].url if ctx.message.attachments else None
        action = action.lower() if action else None
        # --- ADD TAG ---
        if action == 'add':
            resolved_tag_type = 'loop' if (tag_type and tag_type.lower() == 'loop') else 'default'
            if not name:
                return await ctx.send(f'Usage: \"{self.bot.command_prefix}tag add <name> \"content\" [loop]`')
            try:
                await self.tag_manager.add_tag(
                    name=name,
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id,
                    content=content,
                    attachment_url=attachment_url,
                    tag_type=resolved_tag_type
                )
                await ctx.send(f'Tag \"{name}\" (type: {resolved_tag_type}) added successfully.')
            except ValueError as ve:
                await ctx.send(str(ve))
            except Exception as e:
                logger.error(f'Error adding tag: {e}')
                await ctx.send('An error occurred while adding the tag.')
        # --- BORROW TAG ---
        if action == "borrow":
            """
            Usage:
                !tag borrow <tag_name>
                Optionally, specify the original owner: !tag borrow <tag_name> @UserName
            """
            if not name:
                return await ctx.send(
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
                    await ctx.send(
                        f'You have successfully borrowed the tag "{name}" from {owner_display}.'
                    )
                else:
                    await ctx.send(
                        f'You have successfully borrowed the tag "{name}".'
                    )
            except ValueError as ve:
                await ctx.send(str(ve))
            except RuntimeError as re:
                await ctx.send(str(re))
            except Exception as e:
                logger.error(f"Unexpected error during tag borrowing: {e}")
                await ctx.send(
                    "An unexpected error occurred while borrowing the tag."
                )
        # --- LIST TAGS ---
        elif action == 'list':
            filter_tag_type = name.lower() if name and name.lower() in ('loop', 'default') else None
            try:
                tags = await self.tag_manager.list_tags(
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id,
                    tag_type=filter_tag_type
                )
                if not tags:
                    await ctx.send('No tags found.')
                else:
                    tag_list = '\n'.join(f'**{t["name"]}**' for t in tags)
                    # Or, if you still want to display the tag type:
                    # tag_list = '\n'.join(
                    #     f'**{t["name"]}** (type: {t["tag_type"]})'
                    #     for t in tags
                    # )
                    await ctx.send(f'Tags:\n{tag_list}')
            except Exception as e:
                logger.error(f'Error listing tags: {e}')
                await ctx.send('An error occurred while listing your tags.')
        elif action == 'loop':
            if not name:
                return await ctx.send('Usage: \"{self.bot.command_prefix}tag loop on <#channel>` or \"{self.bot.command_prefix}tag loop off`')
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
                    await ctx.send(f'Looping enabled in {channel.mention}.')
                except Exception as e:
                    logger.error(f'Error enabling loop: {e}')
                    await ctx.send('Could not enable loop.')
            elif name.lower() == 'off':
                try:
                    await self.tag_manager.set_loop_config(ctx.guild.id, None, False)
                    self.stop_loop_task()
                    await ctx.send('Looping disabled.')
                except Exception as e:
                    logger.error(f'Error disabling loop: {e}')
                    await ctx.send('Could not disable loop.')
        # --- REMOVE TAG ---
        elif action == 'remove':
            if not name:
                return await ctx.send(f'Usage: \"{self.bot.command_prefix}tag remove <name>`')
            try:
                result = await self.tag_manager.delete_tag(
                    name=name,
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id
                )
                if result > 0:
                    await ctx.send(f'Tag \"{name}\" removed.')
                else:
                    await ctx.send(f'Tag \"{name}\" not found or you do not own it.')
            except Exception as e:
                logger.error(f'Error removing tag: {e}')
                await ctx.send('An error occurred while removing the tag.')
        # --- RENAME TAG ---
        elif action == "rename":
            if not name or not content:
                return await ctx.send(f'Usage: \"{self.bot.command_prefix}tag rename <old_name> <new_name>`')
            old_name = name
            new_name = content
            try:
                row_count = await self.tag_manager.rename_tag(
                    old_name=old_name,
                    new_name=new_name,
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id
                )
                if row_count > 0:
                    await ctx.send(f'Tag \"{old_name}\" renamed to \"{new_name}\".')
                else:
                    await ctx.send(f'Tag \"{old_name}\" not found or you do not own it.')
            except ValueError as ve:
                await ctx.send(str(ve))
            except Exception as e:
                logger.error(f'Error renaming tag: {e}')
                await ctx.send('An error occurred while renaming the tag.')
        # --- UPDATE TAG ---
        elif action == 'update':
            if not name:
                return await ctx.send(f'Usage: {self.bot.command_prefix}tag update <name> \"new content\" [loop|default]`')
            resolved_tag_type = (
                tag_type.lower() if tag_type and tag_type.lower() in ('default', 'loop') else None
            )
            updates = {}
            if content is not None:
                updates['content'] = content
            if attachment_url is not None:
                updates['attachment_url'] = attachment_url
            if resolved_tag_type is not None:
                updates['tag_type'] = resolved_tag_type
            try:
                result = await self.tag_manager.update_tag(
                    name=name,
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id,
                    updates=updates
                )
                if result > 0:
                    await ctx.send(f'Tag \"{name}\" updated.')
                else:
                    await ctx.send(f'Tag \"{name}\" not found or you do not own it.')
            except Exception as e:
                logger.error(f'Error updating tag: {e}')
                await ctx.send('An error occurred while updating the tag.')
        # --- FETCH A TAG BY NAME ---
        else:
            try:
                tag = await self.tag_manager.get_tag(ctx.guild.id, action)
                if tag:
                    content_value = tag.get('content')
                    attachment_url_value = tag.get('attachment_url')
                    if content_value and attachment_url_value:
                        await ctx.send(content_value)
                        await ctx.send(attachment_url_value)
                    elif content_value:
                        await ctx.send(content_value)
                    elif attachment_url_value:
                        await ctx.send(attachment_url_value)
                    else:
                        await ctx.send(f'Tag \"{action}\" has no content.')
                else:
                    await ctx.send(f'Tag \"{action}\" not found.')
            except Exception as e:
                logger.error(f'Error fetching tag \"{action}\": {e}')
                await ctx.send(f'An error occurred while fetching tag \"{action}\".')

#    @commands.hybrid_command(name='tags', description='Display loop tags for the current location.')
#    async def tags(self, ctx: commands.Context):
#        try:
#            if ctx.interaction:
#                async with ctx.typing():
#                    await ctx.interaction.response.defer(ephemeral=True)
#            if not self.predicator.is_release_mode_func(ctx):
#                return
#            location_id = ctx.guild.id
#            tags = await self.tag_manager.list_tags(location_id, tag_type='loop')
#            if not tags:
#                await ctx.send("No loop tags found.")
#                return
#            embeds = []
#            for tag in tags:
#                embed = discord.Embed(
#                    title=f'Loop Tag: {tag["name"]}',
#                    description=tag.get('content', tag.get('attachment_url', 'No content available.')),
#                    color=discord.Color.blurple()
#                )
#                embeds.append(embed)
#            paginator = Paginator(self.bot, ctx, embeds)
#            await paginator.start()
#        except Exception as e:
#            logger.error(f'Error during tag fetching: {e}')

#    @commands.command(name='mod_usage', description=f'Usage: lwipe <all|bot|commands|text|user>')
#    async def mod_usage(self, ctx, limit: int = 1):
#        client = OpenAIUsageClient(api_key=self.config['api_keys']['OpenAI']['api_key'], organization_id=OPENAI_CHAT_HEADERS['OpenAI-Organization'])
#        moderations_usage = await client.get_moderations_usage(
#            start_time=int(time.time()), limit=50
#        )
#        print("Moderations Usage:")
#        for bucket in moderations_usage.data:
#            print(bucket)

    @commands.command(name='wipe', description=f'Usage: lwipe <all|bot|commands|text|user>')
    @commands.has_permissions(manage_messages=True)
    async def wipe(self, ctx, option: str = None, limit: int = 100):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not self.predicator.is_release_mode_func(ctx):
            return
        if limit <= 0 or limit > 100:
            return await ctx.send('Limit must be between 1 and 100.')
        check_function = None
        if option == 'bot':
            check_function = lambda m: m.author == self.bot.user
        elif option == 'all':
            check_function = lambda m: True  # Allow all messages to be deleted
        elif option == 'user':
            user = ctx.message.mentions[0] if ctx.message.mentions else None
            if user:
                check_function = lambda m: m.author == user
            else:
                return await ctx.send('Please mention a user.')
        elif option == 'commands':
            check_function = lambda m: m.content.startswith(ctx.prefix)
        elif option == 'text':
            await ctx.send('Provide text to delete messages containing it.')
            try:
                msg_text = await self.bot.wait_for('message', timeout=30.0, check=lambda m: m.author == ctx.author)
                check_function = lambda m: msg_text.content in m.content
            except asyncio.TimeoutError:
                return await ctx.send('You took too long to provide text. Cancelling operation.')
        else:
            return await ctx.send('Invalid option.')
        total_deleted = 0
        while total_deleted < limit:
            # Purge messages in smaller chunks to avoid hitting rate limits
            deleted = await ctx.channel.purge(limit=min(limit - total_deleted, 10), check=check_function)
            if not deleted:  # Exit loop if no messages were deleted
                break
            total_deleted += len(deleted)
            await asyncio.sleep(1)
        if total_deleted > 0:
            await ctx.send(f'Deleted {total_deleted} messages.')
        else:
            await ctx.send('No messages matched the criteria.')

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

async def setup(bot: commands.Bot):
    await bot.add_cog(Hybrid(bot))
