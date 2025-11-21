-- ==============================
-- 1. Backup existing tables
-- ==============================
-- Optional: rename old tables as backup
ALTER TABLE users RENAME TO users_old;
ALTER TABLE active_text_mutes RENAME TO active_text_mutes_old;
ALTER TABLE active_caps RENAME TO active_caps_old;
ALTER TABLE command_aliases RENAME TO command_aliases_old;
ALTER TABLE active_bans RENAME TO active_bans_old;
ALTER TABLE active_server_voice_mutes RENAME TO active_server_voice_mutes_old;
ALTER TABLE active_cows RENAME TO active_cows_old;
ALTER TABLE active_flags RENAME TO active_flags_old;
ALTER TABLE active_stages RENAME TO active_stages_old;
ALTER TABLE stage_coordinators RENAME TO stage_coordinators_old;
ALTER TABLE active_voice_mutes RENAME TO active_voice_mutes_old;
ALTER TABLE moderation_logs RENAME TO moderation_logs_old;
ALTER TABLE log_channels RENAME TO log_channels_old;
ALTER TABLE role_permissions RENAME TO role_permissions_old;

-- ==============================
-- 2. Create new tables
-- ==============================
CREATE TABLE users (
    discord_snowflake BIGINT PRIMARY KEY,
    moderator_channel_ids BIGINT[],
    coordinator_channel_ids BIGINT[],
    developer_guild_ids BIGINT[],
    server_mute_guild_ids BIGINT[],
    server_muter_guild_ids BIGINT[],
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    coordinator_room_names TEXT[],
    moderator_room_names TEXT[]
);

CREATE TABLE active_text_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    reason TEXT,
    expires_at TIMESTAMPTZ,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name)
);

CREATE TABLE active_caps (
    guild_id BIGINT NOT NULL,
    channel_id BIGINT DEFAULT -1,
    room_name TEXT DEFAULT '',
    moderation_type TEXT NOT NULL CHECK (moderation_type IN ('ban', 'mute', 'tmute')),
    duration TEXT NOT NULL,
    PRIMARY KEY (guild_id, channel_id, room_name, moderation_type)
);

CREATE TABLE command_aliases (
    guild_id   BIGINT NOT NULL,
    alias_type TEXT NOT NULL CHECK (alias_type IN (
        'cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute', 'role', 'unrole'
    )),
    alias_name TEXT NOT NULL,
    channel_id BIGINT DEFAULT -1,
    role_id    BIGINT,
    room_name  TEXT DEFAULT '',
    PRIMARY KEY (guild_id, alias_type, alias_name)
);

CREATE TABLE active_bans (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name)
);

CREATE TABLE active_server_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    channel_id BIGINT DEFAULT -1,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name)
);

CREATE TABLE active_cows (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

CREATE TABLE active_flags (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

CREATE TABLE active_stages (
    guild_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    initiator_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, channel_id, room_name)
);

CREATE TABLE stage_coordinators (
    guild_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    discord_snowflake BIGINT NOT NULL,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name)
);

CREATE TABLE active_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    target TEXT CHECK (target IN ('room', 'user')) DEFAULT 'user',
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name, target)
);

CREATE TABLE moderation_logs (
    id BIGSERIAL PRIMARY KEY,
    action_type TEXT NOT NULL,
    target_discord_snowflake BIGINT,
    executor_discord_snowflake BIGINT NOT NULL,
    guild_id BIGINT NOT NULL,
    channel_id BIGINT,
    reason TEXT,
    metadata JSONB DEFAULT '{}'::jsonb,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE TABLE log_channels (
    guild_id   BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    type       TEXT DEFAULT 'general',
    snowflakes BIGINT[],
    enabled BOOLEAN DEFAULT FALSE,
    PRIMARY KEY (guild_id, channel_id)
);

CREATE TABLE role_permissions (
    role_id BIGINT PRIMARY KEY,
    is_team_member BOOLEAN DEFAULT FALSE
);

CREATE TABLE temporary_rooms (
    guild_snowflake BIGINT NOT NULL,
    room_name TEXT NOT NULL,
    owner_snowflake BIGINT NOT NULL,
    PRIMARY KEY (guild_snowflake, room_name)
);

CREATE TABLE temporary_room_aliases (
    guild_snowflake BIGINT NOT NULL,
    room_name TEXT NOT NULL,
    alias_type TEXT NOT NULL,
    alias_name TEXT NOT NULL,
    role_snowflake BIGINT,
    PRIMARY KEY (guild_snowflake, room_name, alias_type, alias_name)
);

-- ==============================
-- 3. Migrate data from old tables
-- ==============================
INSERT INTO users (discord_snowflake, moderator_channel_ids, coordinator_channel_ids, developer_guild_ids, server_mute_guild_ids, server_muter_guild_ids, updated_at, created_at)
SELECT discord_snowflake, moderator_channel_ids, coordinator_channel_ids, developer_guild_ids, server_mute_guild_ids, server_muter_guild_ids, updated_at, created_at
FROM users_old;

INSERT INTO active_text_mutes (guild_id, discord_snowflake, channel_id, reason, expires_at)
SELECT guild_id, discord_snowflake, channel_id, reason, expires_at
FROM active_text_mutes_old;

INSERT INTO active_caps (guild_id, channel_id, moderation_type, duration)
SELECT guild_id, channel_id, moderation_type, duration
FROM active_caps_old;

INSERT INTO command_aliases (guild_id, alias_type, alias_name, channel_id, role_id)
SELECT guild_id, alias_type, alias_name, channel_id, role_id
FROM command_aliases_old;

INSERT INTO active_bans (guild_id, discord_snowflake, channel_id, expires_at, reason)
SELECT guild_id, discord_snowflake, channel_id, expires_at, reason
FROM active_bans_old;

INSERT INTO active_server_voice_mutes (guild_id, discord_snowflake, expires_at, reason)
SELECT guild_id, discord_snowflake, expires_at, reason
FROM active_server_voice_mutes_old;

INSERT INTO active_cows (guild_id, discord_snowflake, channel_id, created_at)
SELECT guild_id, discord_snowflake, channel_id, created_at
FROM active_cows_old;

INSERT INTO active_flags (guild_id, discord_snowflake, channel_id, expires_at, reason)
SELECT guild_id, discord_snowflake, channel_id, expires_at, reason
FROM active_flags_old;

INSERT INTO active_stages (guild_id, channel_id, initiator_id, expires_at)
SELECT guild_id, channel_id, initiator_id, expires_at
FROM active_stages_old;

INSERT INTO stage_coordinators (guild_id, channel_id, discord_snowflake)
SELECT guild_id, channel_id, discord_snowflake
FROM stage_coordinators_old;

INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason, target)
SELECT guild_id, discord_snowflake, channel_id, expires_at, reason, target
FROM active_voice_mutes_old;

INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason, metadata, created_at)
SELECT action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason, metadata, created_at
FROM moderation_logs_old;

INSERT INTO log_channels (guild_id, channel_id, type, snowflakes, enabled)
SELECT guild_id, channel_id, type, snowflakes, enabled
FROM log_channels_old;

INSERT INTO role_permissions (role_id, is_team_member)
SELECT role_id, is_team_member
FROM role_permissions_old;

-- ==============================
-- 4. Drop old tables if migration is successful
-- ==============================
-- DROP TABLE users_old;
-- DROP TABLE active_text_mutes_old;
-- DROP TABLE active_caps_old;
-- DROP TABLE command_aliases_old;
-- DROP TABLE active_bans_old;
-- DROP TABLE active_server_voice_mutes_old;
-- DROP TABLE active_cows_old;
-- DROP TABLE active_flags_old;
-- DROP TABLE active_stages_old;
-- DROP TABLE stage_coordinators_old;
-- DROP TABLE active_voice_mutes_old;
-- DROP TABLE moderation_logs_old;
-- DROP TABLE log_channels_old;
-- DROP TABLE role_permissions_old;
