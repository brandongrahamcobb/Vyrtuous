import pytest

from vyrtuous.alias.alias_context import AliasContext
from vyrtuous.tests.integration.test_suite import build_message, setup


def test_setup(bot):
    objects = setup(bot)
    msg = build_message(
        author=objects.get("author", None),
        channel=objects.get("text_channel", None),
        content="!aliasname",
        guild=objects.get("guild", None),
        state=objects.get("state", None),
    )
    ctx = AliasContext(message=msg)
    await ctx.setup()
