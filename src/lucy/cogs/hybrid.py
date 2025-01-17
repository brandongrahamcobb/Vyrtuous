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
from lucy.utils.citation_manager import CitationManager
from lucy.utils.combine import combine
from lucy.utils.create_completion import create_completion
from lucy.utils.draw_fingerprint import draw_fingerprint
from lucy.utils.draw_watermarked_molecule import draw_watermarked_molecule
from lucy.utils.game import Game
from lucy.utils.get_mol import get_mol
from lucy.utils.google import google
from lucy.utils.gsrs import gsrs
from lucy.utils.helpers import *
from lucy.utils.paginator import Paginator
from lucy.utils.predicator import Predicator
from lucy.utils.reference_manager import ReferenceManager
from lucy.utils.pdf_manager import PDFManager
from lucy.utils.script import script
from lucy.utils.tag import TagManager
from lucy.utils.unique_pairs import unique_pairs
from lucy.utils.usage import OpenAIUsageClient
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
import pytz
import shlex
import time
import traceback

class Hybrid(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.citation_manager = CitationManager()
        self.config = bot.config
        self.pdf_manager = PDFManager()
        self.predicator = Predicator(self.bot)
        self.reference_manager = ReferencecManager()
        self.tag_manager = TagManager(self.bot.db_pool)
        self.messages = []
        self.game = Game(self.bot, self.bot.db_pool)

    def get_language_code(self, language_name):
        language_name = language_name.lower()
        for lang_code, lang_name in LANGUAGES.items():
            if lang_name.lower() == language_name:
                return lang_code
        return None

    @commands.hybrid_command(name='colorize', description=f'Usage: between `lcolorize 0 0 0` and `lcolorize 255 255 255` or `l colorize <color>`')
    @commands.has_permissions(manage_roles=True) # Do you have manage_roles permissions?
    async def colorize(self, ctx: commands.Context, r: str = commands.parameter(default='blurple', description='Anything between 0 and 255 or a color.'), *, g: str = commands.parameter(default='147', description='Anything betwen 0 and 255.'), b: str = commands.parameter(default='165', description='Anything between 0 and 255.')):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not await self.predicator.is_at_home_func(ctx.guild.id):
            return
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
        quantity: int = commands.parameter(default=1, description='Quantity of glows')
    ):
        try:
            if ctx.interaction:
                async with ctx.typing():
                    await ctx.interaction.response.defer(ephemeral=True)
            if not await self.predicator.is_at_home_func(ctx.guild.id):
                return
            if not self.predicator.is_release_mode_func(ctx):
                return
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

  @commands.hybrid_command(name="addreference", description="Add a new reference.")
    @app_commands.describe(
        location_id="Location ID where the reference is stored.",
        title="Title of the reference.",
        authors="Authors of the reference, separated by commas.",
        publication_year="Publication year of the reference (optional).",
        doi="DOI of the reference (optional).",
        abstract="Abstract of the reference (optional).",
        tags="Tags associated with the reference, separated by commas (optional)."
    )
    async def add_reference(self, ctx: commands.Context, location_id: int, title: str, authors: str, publication_year: Optional[int] = None, doi: Optional[str] = None, abstract: Optional[str] = None, tags: Optional[str] = None):
        user_id = ctx.author.id
        authors_list = [author.strip() for author in authors.split(",")] if authors else []
        tags_list = [int(tag.strip()) for tag in tags.split(",")] if tags else None

        try:
            ref_id = await self.reference_manager.add_reference(
                user_id=user_id,
                location_id=location_id,
                title=title,
                authors=authors_list,
                publication_year=publication_year,
                doi=doi,
                abstract=abstract,
                tags=tags_list
            )
            response = f"✅ Reference added with ID: `{ref_id}`."
            await self.send_response(ctx, response)
        except Exception as e:
            logger.error(f"Error adding reference: {e}")
            await self.send_response(ctx, "❌ Failed to add reference.")

    @commands.hybrid_command(name="deletereference", description="Delete an existing reference.")
    @app_commands.describe(reference_id="The ID of the reference to delete.")
    async def delete_reference(self, ctx: commands.Context, reference_id: int):
        user_id = ctx.author.id

        try:
            success = await self.reference_manager.delete_reference(reference_id, user_id)
            if success:
                response = f"✅ Reference with ID `{reference_id}` deleted."
            else:
                response = f"❌ Reference with ID `{reference_id}` not found or you don't have permission to delete it."
            await self.send_response(ctx, response)
        except Exception as e:
            logger.error(f"Error deleting reference: {e}")
            await self.send_response(ctx, "❌ Failed to delete reference.")

    @commands.hybrid_command(name="listreferences", description="List all your references.")
    @app_commands.describe(
        location_id="Filter by Location ID.",
        tags="Filter by Tag IDs, separated by commas (optional)."
    )
    async def list_references(self, ctx: commands.Context, location_id: int, tags: Optional[str] = None):
        user_id = ctx.author.id
        tags_list = [int(tag.strip()) for tag in tags.split(",")] if tags else None

        try:
            references = await self.reference_manager.list_references(user_id, location_id, tags_list)
            if not references:
                response = "📭 You have no references matching the criteria."
                await self.send_response(ctx, response)
                return

            embed = discord.Embed(title="📑 Your References", color=discord.Color.blue())
            for ref in references:
                authors_str = ', '.join(ref['authors']) if ref['authors'] else 'N/A'
                publication_year = ref.get('publication_year', 'N/A')
                doi = ref.get('doi', 'N/A')
                tags = ref.get('tags', [])
                tags_str = ', '.join(map(str, tags)) if tags else 'N/A'

                embed.add_field(
                    name=f"ID `{ref['id']}`: {ref['title']}",
                    value=(
                        f"**Authors:** {authors_str}\n"
                        f"**Publication Year:** {publication_year}\n"
                        f"**DOI:** {doi}\n"
                        f"**Tags:** {tags_str}"
                    ),
                    inline=False
                )
            await self.send_embed_response(ctx, embed)
        except Exception as e:
            logger.error(f"Error listing references: {e}")
            await self.send_response(ctx, "❌ Failed to list references.")

    @commands.hybrid_command(name="searchreference", description="Search your references by title or authors.")
    @app_commands.describe(
        location_id="Location ID to search within.",
        query_text="Search term for title or authors."
    )
    async def search_reference(self, ctx: commands.Context, location_id: int, query_text: str):
        user_id = ctx.author.id

        try:
            references = await self.reference_manager.search_references(user_id, location_id, query_text)
            if not references:
                response = "🔍 No matching references found."
                await self.send_response(ctx, response)
                return

            embed = discord.Embed(title="🔎 Search Results", color=discord.Color.green())
            for ref in references:
                authors_str = ', '.join(ref['authors']) if ref['authors'] else 'N/A'
                publication_year = ref.get('publication_year', 'N/A')
                doi = ref.get('doi', 'N/A')
                tags = ref.get('tags', [])
                tags_str = ', '.join(map(str, tags)) if tags else 'N/A'

                embed.add_field(
                    name=f"ID `{ref['id']}`: {ref['title']}",
                    value=(
                        f"**Authors:** {authors_str}\n"
                        f"**Publication Year:** {publication_year}\n"
                        f"**DOI:** {doi}\n"
                        f"**Tags:** {tags_str}"
                    ),
                    inline=False
                )
            await self.send_embed_response(ctx, embed)
        except Exception as e:
            logger.error(f"Error searching references: {e}")
            await self.send_response(ctx, "❌ Failed to search references.")

 # ------------- PDF Commands -------------

    @commands.hybrid_command(name="uploadpdf", description="Upload a PDF associated with a reference.")
    @app_commands.describe(
        reference_id="ID of the reference to associate the PDF with.",
        file="The PDF file to upload."
    )
    async def upload_pdf(self, ctx: commands.Context, reference_id: int, file: discord.Attachment):
        user_id = ctx.author.id

        await ctx.defer(ephemeral=True)

        # Check if the user owns the reference
        try:
            reference = await self.reference_manager.get_reference(reference_id, user_id)
            if not reference:
                await self.send_response(ctx, "❌ Reference not found or you don't have permission to modify it.")
                return
        except Exception as e:
            logger.error(f"Error fetching reference: {e}")
            await self.send_response(ctx, "❌ Failed to verify reference.")
            return

        # Ensure the file is a PDF
        if not file.filename.lower().endswith(".pdf"):
            await self.send_response(ctx, "❌ Only PDF files are allowed.")
            return

        try:
            # Download the file
            buffer = io.BytesIO()
            await file.save(buffer)
            buffer.seek(0)

            # Optionally, upload to persistent storage
            # For simplicity, we'll use the Discord CDN URL
            file_url = file.url

            pdf_id = await self.pdf_manager.upload_pdf(reference_id, file_url)
            await self.send_response(ctx, f"✅ PDF uploaded with ID: `{pdf_id}`.")
        except Exception as e:
            logger.error(f"Error uploading PDF: {e}")
            await self.send_response(ctx, "❌ Failed to upload PDF.")

    @commands.hybrid_command(name="annotatepdf", description="Add an annotation to a PDF.")
    @app_commands.describe(
        pdf_id="ID of the PDF to annotate.",
        content="Content of the annotation.",
        page_number="Page number to annotate (optional).",
        highlighted_text="Text to highlight in the PDF (optional)."
    )
    async def annotate_pdf(self, ctx: commands.Context, pdf_id: int, content: str, page_number: Optional[int] = None, highlighted_text: Optional[str] = None):
        user_id = ctx.author.id

        await ctx.defer(ephemeral=True)

        try:
            # Verify if the PDF exists
            pdf = await self.pdf_manager.view_pdf(pdf_id)
            if not pdf:
                await self.send_response(ctx, "❌ PDF not found.")
                return

            # Optionally, verify if the user has access to annotate
            # This depends on your application's logic

            annotation_id = await self.pdf_manager.annotate_pdf(
                pdf_id=pdf_id,
                user_id=user_id,
                page_number=page_number,
                content=content,
                highlighted_text=highlighted_text
            )
            await self.send_response(ctx, f"✅ Annotation added with ID: `{annotation_id}`.")
        except Exception as e:
            logger.error(f"Error annotating PDF: {e}")
            await self.send_response(ctx, "❌ Failed to add annotation.")

    @commands.hybrid_command(name="viewpdf", description="View information about a PDF.")
    @app_commands.describe(pdf_id="ID of the PDF to view.")
    async def view_pdf(self, ctx: commands.Context, pdf_id: int):
        await ctx.defer(ephemeral=True)

        try:
            pdf = await self.pdf_manager.view_pdf(pdf_id)
            if not pdf:
                await self.send_response(ctx, "❌ PDF not found.")
                return

            annotations = await self.pdf_manager.get_annotations(pdf_id)
            annotation_count = len(annotations)

            embed = discord.Embed(title=f"📄 PDF ID: {pdf_id}", color=discord.Color.gold())
            embed.add_field(name="📂 File URL", value=f"[Download PDF]({pdf['file_url']})", inline=False)
            embed.add_field(name="🔗 Reference ID", value=str(pdf['reference_id']), inline=True)
            embed.add_field(name="📝 Annotations", value=f"{annotation_count} found.", inline=True)

            await self.send_embed_response(ctx, embed)
        except Exception as e:
            logger.error(f"Error viewing PDF: {e}")
            await self.send_response(ctx, "❌ Failed to retrieve PDF information.")


    # ------------- Citation Commands -------------

    @commands.hybrid_command(name="generatecitation", description="Generate a citation for a reference.")
    @app_commands.describe(
        reference_id="ID of the reference to cite.",
        citation_style="Citation style (e.g., APA, MLA, Chicago)."
    )
    async def generate_citation(self, ctx: commands.Context, reference_id: int, citation_style: str):
        user_id = ctx.author.id

        await ctx.defer(ephemeral=True)

        try:
            # Verify that the user owns the reference
            reference = await self.reference_manager.get_reference(reference_id, user_id)
            if not reference:
                await self.send_response(ctx, "❌ Reference not found or you don't have permission to cite it.")
                return

            # Call CitationManager to generate the citation
            citation_id = await self.citation_manager.generate_citation(
                reference_id=reference_id,
                user_id=user_id,
                citation_style=citation_style
            )
            await self.send_response(ctx, f"✅ Citation generated with ID: `{citation_id}`.")
        except Exception as e:
            logger.error(f"Error generating citation: {e}")
            await self.send_response(ctx, "❌ Failed to generate citation.")

    @commands.hybrid_command(name="exportbibliography", description="Export your bibliography in a specified citation style.")
    @app_commands.describe(
        citation_style="Citation style for the bibliography (e.g., APA, MLA).",
        location_id="Location ID to export bibliography from."
    )
    async def export_bibliography(self, ctx: commands.Context, citation_style: str, location_id: int):
        user_id = ctx.author.id

        await ctx.defer(ephemeral=True)

        try:
            bibliography = await self.citation_manager.export_bibliography(
                user_id=user_id,
                location_id=location_id,
                citation_style=citation_style
            )

            if not bibliography:
                await self.send_response(ctx, "📄 No citations found to export.")
                return

            # Discord message limit is 2000 characters; if exceeded, send as a file
            if len(bibliography) > 2000:
                buffer = io.StringIO(bibliography)
                file = discord.File(fp=buffer, filename="bibliography.txt")
                await self.send_response(ctx, "📄 Here's your bibliography:", file=file)
            else:
                await self.send_response(ctx, f"📄 **Bibliography ({citation_style}):**\n```\n{bibliography)}\n")

        except Exception as e:
            logger.error(f"Error exporting bibliography: {e}")
            await self.send_response(ctx, "❌ Failed to export bibliography.")

    # ------------- Helper Methods -------------

    async def send_response(self, ctx: commands.Context, content: str, file: discord.File = None):
        """Sends a response to the context, handling both slash and text commands."""
        if isinstance(ctx, commands.Context):
            if ctx.interaction:
                await ctx.interaction.followup.send(content, file=file, ephemeral=True)
            else:
                await ctx.send(content, file=file)
        else:
            await ctx.send(content, file=file)

    async def send_embed_response(self, ctx: commands.Context, embed: discord.Embed):
        """Sends an embed as a response, handling both slash and text commands."""
        if isinstance(ctx, commands.Context):
            if ctx.interaction:
                await ctx.interaction.followup.send(embed=embed, ephemeral=True)
            else:
                await ctx.send(embed=embed)
        else:
            await ctx.send(embed=embed)

    @commands.hybrid_command(name='frame', description='Sends a frame from a number of animal cruelty footage sources.')
    async def frame(self, ctx: commands.Context):
        if ctx.interaction:
            async with ctx.typing():
                await ctx.interaction.response.defer(ephemeral=True)
        if not await self.predicator.is_at_home_func(ctx.guild.id):
            return
        if not self.predicator.is_release_mode_func(ctx):
            return
        video_path = 'frogs.mov'
        output_dir = 'frames'
        frames = extract_random_frames(video_path, output_dir)
        for frame in frames:
            await ctx.send(file=discord.File(frame))

    @commands.command()
    async def leaderboard(self, ctx):
        async with self.db_pool.acquire() as conn:
            rows = await conn.fetch(
                """
                SELECT name, level, exp 
                FROM users 
                ORDER BY level DESC, exp DESC 
                LIMIT 10
                """
            )
        if not rows:
            await ctx.send("The leaderboard is empty! Be the first to interact with the bot and gain XP!")
            return
        leaderboard_message = "**🏆 Leaderboard 🏆**\n"
        for i, row in enumerate(rows, start=1):
            leaderboard_message += f"{i}. **{row['name']}** - Level {row['level']}, {row['exp']:.2f} XP\n"
        await ctx.send(leaderboard_message)

    @commands.command()
    async def level(self, ctx, member: discord.Member = None):
        member = member or ctx.author
        user_id = member.id
        if user_id not in self.game.users:
            await ctx.send(f"{member.mention} has not interacted with the bot yet.")
            return
        user = self.game.users[user_id]
        xp_needed = self.game.get_xp_for_level(user["level"] + 1) - user["xp"]
        await ctx.send(f"{member.mention} is at level {user['level']} with {user['xp']:.2f} XP. "
                       f"You need {xp_needed:.2f} XP to reach level {user['level'] + 1}.")

#    @commands.command()
#    async def languages(self, ctx):
#        supported_languages = ', '.join(LANGUAGES.values())
#        await ctx.send(f'Supported languages are:\n{supported_languages}')
#
    @commands.command(name='script', description=f'Usage: lscript <NIV/ESV/Quran> <Book>.<Chapter>.<Verse>')
    async def script(self, ctx: commands.Context, version: str, *, reference: str):
        try:
            if ctx.interaction:
                async with ctx.typing():
                    await ctx.interaction.response.defer(ephemeral=True)
            if not await self.predicator.is_at_home_func(ctx.guild.id):
                return
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
        if not await self.predicator.is_at_home_func(ctx.guild.id):
            return
        if not self.predicator.is_release_mode_func(ctx):
            return
        results = google(query)
        embed = discord.Embed(title=f'Search Results for \"{query}\"', color=discord.Color.blue())
        for result in results:
            embed.add_field(name=result['title'], value=result['link'], inline=False)
        await ctx.send(embed=embed)

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
                    # Change here: only join the tag's name (plus type if desired)
                    # For name ONLY:
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

    @commands.hybrid_command(name='tags', description='Display loop tags for the current location.')
    async def tags(self, ctx: commands.Context):
        try:
            if ctx.interaction:
                async with ctx.typing():
                    await ctx.interaction.response.defer(ephemeral=True)
            if not self.predicator.is_release_mode_func(ctx):
                return
            location_id = ctx.guild.id
            tags = await self.tag_manager.list_tags(location_id, tag_type='loop')
            if not tags:
                await ctx.send("No loop tags found.")
                return
            embeds = []
            for tag in tags:
                embed = discord.Embed(
                    title=f'Loop Tag: {tag["name"]}',
                    description=tag.get('content', tag.get('attachment_url', 'No content available.')),
                    color=discord.Color.blurple()
                )
                embeds.append(embed)
            paginator = Paginator(self.bot, ctx, embeds)
            await paginator.start()
        except Exception as e:
            logger.error(f'Error during tag fetching: {e}')

    @commands.command(name='mod_usage', description=f'Usage: lwipe <all|bot|commands|text|user>')
    async def mod_usage(self, ctx, limit: int = 1):
        client = OpenAIUsageClient(api_key=self.config['api_keys']['OpenAI']['api_key'], organization_id=OPENAI_CHAT_HEADERS['OpenAI-Organization'])
        moderations_usage = await client.get_moderations_usage(
            start_time=int(time.time()), limit=50
        )
        print("Moderations Usage:")
        for bucket in moderations_usage.data:
            print(bucket)

    @commands.command(name='wipe', description=f'Usage: lwipe <all|bot|commands|text|user>')
    @commands.has_permissions(manage_messages=True)
    async def wipe(self, ctx, option: str = None, limit: int = 100):
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
