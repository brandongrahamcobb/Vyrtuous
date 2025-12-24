
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

CREATE TABLE administrators (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    role_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, member_snowflake, role_snowflake)
);

DROP TABLE active_server_voice_mutes;
CREATE TABLE active_server_voice_mutes (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, member_snowflake)
);

ALTER TABLE active_caps RENAME TO active_caps_old;
CREATE TABLE active_caps (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    duration_seconds INTEGER NOT NULL,
    guild_snowflake BIGINT NOT NULL,
    channel_snowflake BIGINT DEFAULT -1,
    moderation_type TEXT NOT NULL CHECK (moderation_type IN ('ban', 'mute', 'tmute')),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, channel_snowflake, moderation_type)
);
INSERT INTO active_caps (duration_seconds, guild_snowflake, channel_snowflake, moderation_type)
SELECT duration_seconds, guild_id, channel_id, moderation_type
FROM active_caps_old;
DROP TABLE active_caps_old;

ALTER TABLE active_text_mutes RENAME TO active_text_mutes_old;
CREATE TABLE active_text_mutes (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_at TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);
INSERT INTO active_text_mutes (channel_snowflake, guild_snowflake, member_snowflake, reason, expires_at)
SELECT channel_id, guild_id, discord_snowflake, reason, expires_at
FROM active_text_mutes_old;
DROP TABLE active_text_mutes_old;

ALTER TABLE active_bans RENAME TO active_bans_old;
CREATE TABLE active_bans (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_at TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);
INSERT INTO active_bans (channel_snowflake, guild_snowflake, member_snowflake, reason, expires_at)
SELECT channel_id, guild_id, discord_snowflake, reason, expires_at
FROM active_bans_old;
DROP TABLE active_bans_old;

ALTER TABLE active_voice_mutes RENAME TO active_voice_mutes_old;
CREATE TABLE active_voice_mutes (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_at TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    target TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);
INSERT INTO active_voice_mutes (channel_snowflake, guild_snowflake, member_snowflake, reason, expires_at)
SELECT channel_id, guild_id, discord_snowflake, reason, expires_at
FROM active_voice_mutes_old;
DROP TABLE active_voice_mutes_old;

ALTER TABLE temporary_rooms RENAME TO temporary_rooms_old;
CREATE TABLE temporary_rooms (
    channel_snowflake BIGINT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    room_name TEXT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, room_name)
);
INSERT INTO temporary_rooms (
    channel_snowflake,
    guild_snowflake,
    member_snowflake,
    room_name
)
SELECT
    room_snowflake,
    guild_snowflake,
    owner_snowflake,
    room_name
FROM temporary_rooms_old;
DROP TABLE temporary_rooms_old;

DROP TABLE active_stages;
CREATE TABLE active_stages (
    channel_snowflake BIGINT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_at TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    PRIMARY KEY (guild_snowflake, channel_snowflake)
);
