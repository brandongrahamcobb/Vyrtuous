def test_read():
    read_the_docs()


def read_the_docs():
    with open("tests/unit_tests/docs.python.org", "r") as f:
        print(f.read())


if __name__ == "__main__":
    test_read()

    @commands.command(name="overwrites")
    async def overwrites(self, ctx, channel: discord.abc.GuildChannel | None = None):
        channel = channel or ctx.channel
        member_count = 0
        role_count = 0
        total_count = 0
        for target, overwrite in channel.overwrites.items():
            if any(v is not None for v in overwrite):
                total_count += 1
                if isinstance(target, discord.Member):
                    member_count += 1
                else:
                    role_count += 1
        await ctx.send(
            f"Channel: {channel.name}\nTotal overwrites: {total_count}\nRole overwrites: {role_count}\nMember overwrites: {member_count}"
        )
