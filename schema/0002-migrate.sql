
CREATE TABLE moderators (
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    channel_snowflakes BIGINT NOT NULL
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake
);

INSERT INTO moderators (channel_snowflake, guild_snowflake, member_snowflake)
SELECT
    channel_id,
    801609515391778826,
    discord_snowflake
FROM users
CROSS JOIN unnest(COALESCE(moderator_channel_ids, '{}')) AS channel_id;

CREATE TABLE coordinators (
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    channel_snowflakes BIGINT NOT NULL
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake
);

INSERT INTO coordinators (channel_snowflake, guild_snowflake, member_snowflake)
SELECT
    channel_id,
    801609515391778826,
    discord_snowflake
FROM users
CROSS JOIN unnest(COALESCE(coordinator_channel_ids, '{}')) AS channel_id;
