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
