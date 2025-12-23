
CREATE TABLE moderators (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

INSERT INTO moderators (channel_snowflake, guild_snowflake, member_snowflake)
SELECT
    channel_id,
    801609515391778826,
    discord_snowflake
FROM users
CROSS JOIN unnest(COALESCE(moderator_channel_ids, '{}')) AS channel_id;

CREATE TABLE coordinators (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

INSERT INTO coordinators (channel_snowflake, guild_snowflake, member_snowflake)
SELECT
    channel_id,
    801609515391778826,
    discord_snowflake
FROM users
CROSS JOIN unnest(COALESCE(coordinator_channel_ids, '{}')) AS channel_id;

CREATE TABLE developers (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, member_snowflake)
);
