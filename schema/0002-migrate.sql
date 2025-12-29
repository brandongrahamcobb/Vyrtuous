
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
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, channel_snowflake)
);

ALTER TABLE command_aliases RENAME TO command_aliases_old;
CREATE TABLE command_aliases (
    alias_name TEXT NOT NULL,
    alias_type TEXT NOT NULL CHECK (alias_type IN (
        'cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute', 'role', 'unrole'
    )),
    channel_snowflake BIGINT DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    role_snowflake BIGINT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (alias_name, alias_type, guild_snowflake)
);
DROP TABLE command_aliases_old;

INSERT INTO command_aliases (
    alias_name,
    alias_type,
    channel_snowflake,
    guild_snowflake,
    role_snowflake
)
SELECT
    alias_name,
    alias_type,
    channel_id,
    guild_id,
    role_id
FROM command_aliases_old;
DROP TABLE command_aliases_old;

ALTER TABLE statistic_channels RENAME TO statistic_channels_old;

CREATE TABLE statistic_channels (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    enabled BOOLEAN DEFAULT FALSE,
    guild_snowflake BIGINT NOT NULL,
    snowflakes BIGINT[],
    statistic_type TEXT DEFAULT 'general',
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake)
);
INSERT INTO statistic_channels (
    channel_snowflake,
    enabled,
    guild_snowflake,
    snowflakes,
    statistic_type
)
SELECT
    channel_id,
    enabled,
    guild_id,
    snowflakes,
    type
FROM statistic_channels_old;

BEGIN;

-- 1. Rename the old table
ALTER TABLE active_cows RENAME TO active_cows_old;

-- 2. Create the new table
CREATE TABLE vegans (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

-- 3. Copy data from the old table
INSERT INTO vegans (channel_snowflake, guild_snowflake, member_snowflake, created_at)
SELECT channel_id, guild_id, discord_snowflake, created_at
FROM active_cows_old;

-- 4. Drop the old table
DROP TABLE active_cows_old;

COMMIT;


ALTER TABLE active_flags RENAME TO active_flags_old;
CREATE TABLE active_flags (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_at TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);
INSERT INTO active_flags (channel_snowflake, guild_snowflake, member_snowflake, reason, expires_at)
SELECT channel_id, guild_id, discord_snowflake, reason, expires_at
FROM active_flags_old;
DROP TABLE active_flags_old;

BEGIN TRANSACTION;

-- 1. Rename the existing table
ALTER TABLE command_aliases RENAME TO command_aliases_old;

-- 2. Create the new table with updated column names and CHECK constraint
CREATE TABLE command_aliases (
    alias_type        TEXT NOT NULL CHECK (alias_type IN (
        'vegan', 'carnist', 'voice_mute', 'unvoice_mute', 'ban', 'unban', 'flag', 'unflag', 'text_mute', 'untext_mute', 'role', 'unrole'
    )),
    alias_name        TEXT NOT NULL,
    channel_snowflake BIGINT DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake   BIGINT NOT NULL,
    role_snowflake    BIGINT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (alias_name, alias_type, guild_snowflake)
);

-- 3. Copy and rename alias_type in existing data, mapping old column names to new ones
INSERT INTO command_aliases (guild_snowflake, alias_type, alias_name, channel_snowflake, role_snowflake)
SELECT
    guild_id,
    CASE alias_type
        WHEN 'mute' THEN 'voice_mute'
        WHEN 'unmute' THEN 'unvoice_mute'
        WHEN 'tmute' THEN 'text_mute'
        WHEN 'untmute' THEN 'untext_mute'
        WHEN 'cow' THEN 'vegan'
        WHEN 'uncow' THEN 'carnist'
        ELSE alias_type
    END AS alias_type,
    alias_name,
    channel_id,
    role_id
FROM command_aliases_old;

-- 4. Drop the old table
DROP TABLE command_aliases_old;

COMMIT;

