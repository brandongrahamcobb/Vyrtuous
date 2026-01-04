BEGIN;
DROP TABLE IF EXISTS statistic_channels_old;
ALTER TABLE statistic_channels RENAME TO statistic_channels_old;
CREATE TABLE IF NOT EXISTS history (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ,
    enabled BOOLEAN DEFAULT FALSE,
    entry_type TEXT DEFAULT 'general',
    guild_snowflake BIGINT NOT NULL,
    snowflakes BIGINT[],
    updated_at TIMESTAMPTZ,
    PRIMARY KEY (channel_snowflake, guild_snowflake)
);
INSERT INTO history (channel_snowflake, created_at, enabled, entry_type, guild_snowflake, snowflakes, updated_at)
SELECT channel_snowflake, created_at, enabled, statistic_type, guild_snowflake, snowflakes, updated_at
FROM statistic_channels_old
ON CONFLICT (channel_snowflake, guild_snowflake) DO UPDATE
SET enabled = EXCLUDED.enabled,
    entry_type = EXCLUDED.entry_type,
    snowflakes = EXCLUDED.snowflakes,
    updated_at = EXCLUDED.updated_at;
COMMIT;


ALTER TABLE moderation_logs RENAME TO moderation_logs_old;

CREATE TABLE moderation_logs (
    action_type TEXT NOT NULL,
    channel_snowflake BIGINT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    executor_member_snowflake BIGINT,
    expires_at TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    highest_role TEXT,
    is_modification BOOLEAN NOT NULL DEFAULT FALSE,
    reason TEXT,
    target_member_snowflake BIGINT,
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

INSERT INTO moderation_logs (
    action_type,
    channel_snowflake,
    created_at,
    executor_member_snowflake,
    expires_at,
    guild_snowflake,
    highest_role,
    is_modification,
    reason,
    target_member_snowflake,
    updated_at
)
SELECT
    action_type,
    channel_id,
    NOW(),
    executor_discord_snowflake,
    NULL,
    guild_id,
    NULL,
    FALSE,
    reason,
    target_discord_snowflake,
    NOW()
FROM moderation_logs_old;

