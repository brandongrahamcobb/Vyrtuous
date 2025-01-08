import asyncio

class Paginator:
    def __init__(self, bot, ctx, pages):
        """
        A custom paginator for Discord embeds.
        :param bot: The bot instance
        :param ctx: The context of the command
        :param pages: A list of discord.Embed objects
        """
        self.bot = bot
        self.ctx = ctx
        self.pages = pages
        self.current_page = 0
        self.message = None

    async def start(self):
        """
        Starts the pagination by sending the first embed and adding reaction controls.
        """
        if not self.pages:
            await self.ctx.send("There are no tags to display.")
            return

        # Send the first page
        self.message = await self.ctx.send(embed=self.pages[self.current_page])

        # Add reactions for navigation
        await self.message.add_reaction("⬅️")
        await self.message.add_reaction("➡️")
        await self.message.add_reaction("⏹️")

        def check(reaction, user):
            return (
                user == self.ctx.author
                and reaction.message.id == self.message.id
                and str(reaction.emoji) in ["⬅️", "➡️", "⏹️"]
            )

        while True:
            try:
                # Wait for a reaction
                reaction, user = await self.bot.wait_for("reaction_add", timeout=60.0, check=check)

                if str(reaction.emoji) == "⬅️":
                    # Go to the previous page
                    if self.current_page > 0:
                        self.current_page -= 1
                        await self.message.edit(embed=self.pages[self.current_page])
                elif str(reaction.emoji) == "➡️":
                    # Go to the next page
                    if self.current_page < len(self.pages) - 1:
                        self.current_page += 1
                        await self.message.edit(embed=self.pages[self.current_page])
                elif str(reaction.emoji) == "⏹️":
                    # Stop pagination
                    await self.message.clear_reactions()
                    break

                # Remove the user's reaction
                await self.message.remove_reaction(reaction.emoji, user)
            except asyncio.TimeoutError:
                # Clear reactions if the user takes too long
                await self.message.clear_reactions()
                break

