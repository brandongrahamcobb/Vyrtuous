from datetime import datetime, timezone
from vyrtuous.inc.helpers import PRIVILEGED_AUTHOR_ID, VOICE_CHANNEL_ONE_ID
from vyrtuous.database.database_factory import DatabaseFactory


set_kwargs = {"created_at": datetime.now(timezone.utc), "reason": "Test reason"}

where_kwargs = {
    "channel_snowflake": VOICE_CHANNEL_ONE_ID,
    "member_snowflake": PRIVILEGED_AUTHOR_ID,
}


def test_update():
    pass
