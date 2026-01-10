import pytest
from contextlib import ExitStack
from unittest.mock import patch

from vyrtuous.utils.search import Search
from vyrtuous.config import Config
from vyrtuous.inc.helpers import MESSAGE_ID, \
    GUILD_ID, \
    PRIVILEGED_AUTHOR_ID, \
    ROLE_ID, \
    VOICE_CHANNEL_ONE_ID
from vyrtuous.service.state_service import State
from vyrtuous.tests.black_box.make_mock_objects import create_message
from vyrtuous.tests.black_box.test_suite import bot, \
    capture, \
    guild, \
    privileged_author, \
    prepare_context, \
    prepare_discord_state, \
    text_channel
from vyrtuous.enhanced_member.administrator import Administrator
from vyrtuous.enhanced_member.coordinator import Coordinator
from vyrtuous.enhanced_member.developer import Developer
from vyrtuous.enhanced_member.moderator import Moderator
from vyrtuous.moderation_action.ban import Ban
from vyrtuous.moderation_action.flag import Flag
from vyrtuous.moderation_action.text_mute import TextMute
from vyrtuous.moderation_action.vegan import Vegan
from vyrtuous.moderation_action.voice_mute import VoiceMute

@pytest.mark.asyncio
async def test_resolve_objs_from_target(bot, guild, privileged_author, text_channel):
   
    config = Config.get_config()

    obj_classes = {
        'enhanced_member': [
            Administrator,
            Coordinator,
            Developer,
            Moderator
        ],
        'moderation_action': [
            Ban,
            Flag,
            TextMute,
            Vegan,
            VoiceMute
        ]
    }

    targets = [
        'all',
        VOICE_CHANNEL_ONE_ID,
        ROLE_ID,
        GUILD_ID,
        PRIVILEGED_AUTHOR_ID
        
    ]

    content = ''
    highest_role = 'System Owner'
    prefix = config['discord_command_prefix']

    message = create_message(
        allowed_mentions=True,
        author=privileged_author,
        bot=bot,
        channel=text_channel,
        content=content,
        guild=guild,
        id=MESSAGE_ID
    )
    mock_bot_user = guild.me
    with patch.object(bot, "_connection", create=True) as mock_conn, ExitStack() as stack, prepare_discord_state(privileged_author, bot, text_channel, content, guild, highest_role):
        mock_conn.user = mock_bot_user
        mock_conn.return_value = mock_bot_user
        ctx = await prepare_context(bot=bot, message=message, prefix=prefix)
        state = State(ctx)
        async with capture(privileged_author, text_channel) as captured:
            for group, classes in obj_classes.items():
                for item in classes:
                    for target in targets:
                        objs = \
                            await Search.resolve_objs_from_target(
                                ctx_interaction_or_message=ctx,
                                obj=item,
                                state=state,
                                target=target
                            )
                        if not isinstance(objs, list):
                            print(objs.content)
                            continue
                        for obj in objs:
                            match target:
                                case 'all':
                                    assert True
                                case _ if target == GUILD_ID:
                                    assert obj.guild_snowflake == VOICE_CHANNEL_ONE_ID
                                case _ if target == ROLE_ID:
                                    if hasattr(obj, 'role_snowflake'):
                                        assert obj.role_snowflake == VOICE_CHANNEL_ONE_ID
                                    else:
                                        assert True
                                case _ if target == VOICE_CHANNEL_ONE_ID:
                                    assert obj.channel_snowflake == VOICE_CHANNEL_ONE_ID

