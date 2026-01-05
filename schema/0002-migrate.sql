DROP TABLE history;
CREATE TABLE history (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    enabled BOOLEAN DEFAULT FALSE,
    guild_snowflake BIGINT NOT NULL,
    id BIGSERIAL PRIMARY KEY,
    snowflakes BIGINT[],
    statistic_type TEXT DEFAULT 'general',
    updated_at TIMESTAMPTZ DEFAULT NOW()
);
