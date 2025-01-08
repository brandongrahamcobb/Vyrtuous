''' hybrid.py The purpose of this program is to be an extension to a Discord bot to provide the command functionality from cd lucy.lucy./lucy.lucy./.
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
from discord import Embed
from discord.ext import commands, tasks
from PIL import Image
from random import randint
from typing import Optional
from lucy.utils.frames import extract_random_frames
from lucy.utils.add_watermark import add_watermark
from lucy.utils.average_score import average_score
from lucy.utils.backup import perform_backup, setup_backup_directory
from lucy.utils.combine import combine
from lucy.utils.create_completion import create_completion
from lucy.utils.draw_fingerprint import draw_fingerprint
from lucy.utils.draw_watermarked_molecule import draw_watermarked_molecule
from lucy.utils.get_mol import get_mol
from lucy.utils.google import google
from lucy.utils.gsrs import gsrs
from lucy.utils.helpers import *
from lucy.utils.script import script
from lucy.utils.unique_pairs import unique_pairs
from lucy.utils.tag import TagManager
from random import choice
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
import traceback

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
        self.loop_task = None
        self.daily_loop.start()  # Start the daily check loop when cog is loaded
        self.channel_guild_map: Dict[int, int] = {
            798967615636504657: 730907954345279591,
            730907954877956179: 730907954345279591,
        }
        self.guild_loops_index = defaultdict(int)

    def cog_unload(self):
        self.daily_loop.cancel()

    @tasks.loop(minutes=1)
    async def daily_loop(self):
        """
        Checks every minute if it's 10:00 PM EST.
        If it is, picks the "next" loop tag for each channel in each guild, in order.
        """
        await self.bot.wait_until_ready()

        # Current time in EST
        est_tz = pytz.timezone('US/Eastern')
        now_est = datetime.datetime.now(est_tz)

        # If it's 10:00 PM (22:00) local EST time, do our sends
        if now_est.hour == 22 and now_est.minute == 0:
            # Group channels by guild
            guild_channels_map = {}
            for channel_id, guild_id in self.channel_guild_map.items():
                guild_channels_map.setdefault(guild_id, []).append(channel_id)

            # For each guild, fetch the list of loop tags.
            for guild_id, channel_ids in guild_channels_map.items():
                loop_tags = await self.tag_manager.list_tags(
                    location_id=guild_id,
                    tag_type='loop'
                )

                # Filter out any tags that have neither content nor an attachment
                loop_tags = [
                    t for t in loop_tags
                    if t.get('content') or t.get('attachment_url')
                ]

                if not loop_tags:
                    # If no loop tags, let each channel know or just skip
                    for cid in channel_ids:
                        channel = self.bot.get_channel(cid)
                        if channel:
                            await channel.send("No loop tags found for this guild.")
                    # Move on to the next guild
                    continue

                # Grab the current index for this guild
                current_index = self.guild_loops_index[guild_id]

                # For each channel (in order), pick the next tag
                for cid in channel_ids:
                    channel = self.bot.get_channel(cid)
                    if channel:
                        tag = loop_tags[current_index % len(loop_tags)]
                        msg = tag['content'] or tag['attachment_url']
                        if msg:
                            await channel.send(msg)
                        # Move to the next index for *each channel*
                        current_index += 1

                # After we've used a few tags for these channels,
                # update the guild's index
                self.guild_loops_index[guild_id] = current_index

    @daily_loop.before_loop
    async def before_daily_loop(self):
        """Just an optional hook before the loop starts."""
        print("Daily loop is waiting until bot is ready...")

    @commands.hybrid_command(
        name='tag',
        description='Manage or retrieve tags. Sub-actions: add, update, remove, list, loop.'
    )
    async def tag_command(
        self,
        ctx: commands.Context,
        action: str = commands.parameter(default=None, description='Action: add, update, remove, list, loop.'),
        name: Optional[str] = commands.parameter(default=None, description='Name of the tag.'),
        content: Optional[str] = commands.parameter(default=None, description='Content for the tag (if applicable).'),
        tag_type: Optional[str] = commands.parameter(default=None, description='Optional tag type: default or loop.')
    ):
        """
        Examples:
        !tag add greet "Hello, world!"
        !tag update greet "Hello again!"
        !tag remove greet
        !tag list
        !tag list loop
        !tag loop on #some-channel
        !tag loop off
        !tag greet
        """
        # If there's an attachment, store the first attachment's URL
        attachment_url = ctx.message.attachments[0].url if ctx.message.attachments else None
        action = action.lower() if action else None

        # --- ADD TAG ---
        if action == "add":
            resolved_tag_type = 'loop' if (tag_type and tag_type.lower() == 'loop') else 'default'
            if not name:
                return await ctx.send("Usage: `!tag add <name> \"content\" [loop]`")
            try:
                await self.tag_manager.add_tag(
                    name=name,
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id,
                    content=content,
                    attachment_url=attachment_url,
                    tag_type=resolved_tag_type
                )
                await ctx.send(f"Tag `{name}` (type: {resolved_tag_type}) added successfully.")
            except ValueError as ve:
                await ctx.send(str(ve))
            except Exception as e:
                logger.error(f"Error adding tag: {e}")
                await ctx.send("An error occurred while adding the tag.")

        # --- UPDATE TAG ---
        elif action == "update":
            if not name:
                return await ctx.send("Usage: `!tag update <name> \"new content\" [loop|default]`")
            resolved_tag_type = (
                tag_type.lower() if tag_type and tag_type.lower() in ('default', 'loop') else None
            )

            # Build a dict of fields to update
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
                    await ctx.send(f"Tag `{name}` updated.")
                else:
                    await ctx.send(f"Tag `{name}` not found or you do not own it.")
            except Exception as e:
                logger.error(f"Error updating tag: {e}")
                await ctx.send("An error occurred while updating the tag.")

        # --- REMOVE TAG ---
        elif action == "remove":
            if not name:
                return await ctx.send("Usage: `!tag remove <name>`")
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

        # --- LIST TAGS ---
        elif action == "list":
            # If user typed `!tag list loop`, filter by loop
            # If user typed `!tag list default`, filter by default
            # Otherwise, list them all
            filter_tag_type = name.lower() if name and name.lower() in ('loop', 'default') else None
            try:
                tags = await self.tag_manager.list_tags(
                    location_id=ctx.guild.id,
                    owner_id=ctx.author.id,
                    tag_type=filter_tag_type
                )
                if not tags:
                    await ctx.send("No tags found.")
                else:
                    tag_list = "\n".join(
                        f"**{t['name']}** (type: {t['tag_type']}): {t['content'] or t['attachment_url']}"
                        for t in tags
                    )
                    await ctx.send(f"Tags:\n{tag_list}")
            except Exception as e:
                logger.error(f"Error listing tags: {e}")
                await ctx.send("An error occurred while listing your tags.")

        # --- LOOP TAGS (on/off) ---
        elif action == "loop":
            if not name:
                return await ctx.send("Usage: `!tag loop on <#channel>` or `!tag loop off`")

            if name.lower() == "on":
                # If user typed: !tag loop on #channel
                # content might contain "<#1234567890>"
                # fallback to current channel if user didn't mention any
                channel = ctx.channel
                if content and content.startswith("<#") and content.endswith(">"):
                    channel_id = int(content.strip("<#>"))
                    maybe_chan = self.bot.get_channel(channel_id)
                    if maybe_chan is not None:
                        channel = maybe_chan

                try:
                    await self.tag_manager.set_loop_config(ctx.guild.id, channel.id, True)
                    self.start_loop_task(channel)
                    await ctx.send(f"Looping enabled in {channel.mention}.")
                except Exception as e:
                    logger.error(f"Error enabling loop: {e}")
                    await ctx.send("Could not enable loop.")

            elif name.lower() == "off":
                try:
                    await self.tag_manager.set_loop_config(ctx.guild.id, None, False)
                    self.stop_loop_task()
                    await ctx.send("Looping disabled.")
                except Exception as e:
                    logger.error(f"Error disabling loop: {e}")
                    await ctx.send("Could not disable loop.")

        # --- FETCH A TAG BY NAME ---
        else:
            # If user does something like !tag greet
            # action is actually the tag name, so we fetch it.
            try:
                tag = await self.tag_manager.get_tag(ctx.guild.id, action)
                if tag:
                    content_value = tag.get("content")
                    attachment_url_value = tag.get("attachment_url")
                    if content_value and attachment_url_value:
                        await ctx.send(content_value)
                        await ctx.send(attachment_url_value)
                    elif content_value:
                        await ctx.send(content_value)
                    elif attachment_url_value:
                        await ctx.send(attachment_url_value)
                    else:
                        await ctx.send(f"Tag `{action}` has no content.")
                else:
                    await ctx.send(f"Tag `{action}` not found.")
            except Exception as e:
                logger.error(f"Error fetching tag '{action}': {e}")
                await ctx.send(f"An error occurred while fetching tag `{action}`.")

    def start_loop_task(self, channel: discord.TextChannel):
        """
        Start a background task that periodically sends random "loop" tags.
        """
        if self.loop_task is None or self.loop_task.done():
            self.loop_task = asyncio.create_task(self.loop_tags(channel))

    def stop_loop_task(self):
        """
        Stop the background loop (if any).
        """
        if self.loop_task and not self.loop_task.done():
            self.loop_task.cancel()
            self.loop_task = None

    async def loop_tags(self, channel: discord.TextChannel):
        """
        Periodically picks random 'loop' tags from the DB for this guild and sends them.
        """
        while True:
            try:
                # Only loop tags of type 'loop' for the guild
                loop_tags = await self.tag_manager.list_tags(channel.guild.id, tag_type='loop')
                if loop_tags:
                    random_tag = choice(loop_tags)
                    message_text = random_tag['content'] or random_tag['attachment_url'] or ''
                    if message_text:
                        await channel.send(message_text)

                # Sleep 5 minutes between looped messages
                await asyncio.sleep(300)
            except asyncio.CancelledError:
                # Task canceled externally (stop_loop_task called)
                break
            except Exception as e:
                logger.error(f"Error during loop_tags: {e}")

    @commands.command(name='script', description='Usage !script <NIV/ESV> <Book>.<Chapter>.<Verse>', hidden=True)
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

    @commands.hybrid_command(name='search', description='Usage: !search <query>. Search Google.', hidden=True)
    @at_home()
    async def search(self, ctx: commands.Context, *, query: str = commands.parameter(default=None, description='Google search a query.')):
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
        results = google(query)
        embed = discord.Embed(title=f'Search Results for `{query}`', color=discord.Color.blue())
        for result in results:
            embed.add_field(name=result['title'], value=result['link'], inline=False)
        await ctx.send(embed=embed)

    @commands.hybrid_command(name='frame', description='', hidden=True)
    @at_home()
    async def frame(self, ctx: commands.Context):
        video_path = 'frogs.mov'
        output_dir = 'frames'
        frames = extract_random_frames(video_path, output_dir)
        for frame in frames:
            await ctx.send(file=discord.File(frame))

async def setup(bot: commands.Bot):
    await bot.add_cog(Hybrid(bot))
