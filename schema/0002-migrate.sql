ALTER TABLE active_bans DROP COLUMN display_name;
ALTER TABLE active_bans DROP COLUMN expired;
ALTER TABLE active_bans ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_caps ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_flags DROP COLUMN display_name;
ALTER TABLE active_flags DROP COLUMN expires_in;
ALTER TABLE active_flags ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_server_voice_mutes DROP COLUMN display_name;
ALTER TABLE active_stages DROP COLUMN expired;
ALTER TABLE active_stages RENAME TO active_automute_channels;
ALTER TABLE active_text_mutes DROP COLUMN expired;
ALTER TABLE active_text_mutes DROP COLUMN role_snowflake;
ALTER TABLE active_text_mutes DROP COLUMN display_name;
ALTER TABLE active_text_mutes ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_voice_mutes DROP COLUMN expired;
ALTER TABLE active_voice_mutes DROP COLUMN display_name;
ALTER TABLE active_voice_mutes ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_voice_mutes
DROP CONSTRAINT active_voice_mutes_pkey;
CREATE UNIQUE INDEX active_voice_mutes_unique_idx
ON active_voice_mutes (
    guild_snowflake,
    member_snowflake,
    target,
    channel_snowflake
)
NULLS NOT DISTINCT;
ALTER TABLE administrators DROP COLUMN display_name;
DROP TABLE ban_roles;
DELETE FROM command_aliases
WHERE category = 'vegan';
ALTER TABLE command_aliases
DROP CONSTRAINT command_aliases_category_check;
ALTER TABLE command_aliases
ADD CONSTRAINT command_aliases_category_check
CHECK (
    category = ANY (
        ARRAY[
            'vmute',
            'ban',
            'flag',
            'tmute',
            'role'
        ]::text[]
    )
);
ALTER TABLE coordinators DROP COLUMN display_name;
ALTER TABLE developers DROP COLUMN display_name;
DROP TABLE guild_owners;
DROP TABLE hide_roles;
ALTER TABLE moderators DROP COLUMN display_name;
ALTER TABLE moderation_logs ADD COLUMN target TEXT DEFAULT 'user';
ALTER TABLE moderation_logs ADD COLUMN role_snowflake BIGINT;
DROP TABLE sysadmin;
DROP TABLE temporary_blacklist;
DROP TABLE text_mute_roles;
DROP TABLE users;
ALTER TABLE vegans DROP COLUMN display_name;
DROP TRIGGER set_expired_active_bans ON active_bans;
DROP TRIGGER set_expired_active_stages ON active_automute_channels;
DROP TRIGGER set_expired_active_text_mutes ON active_text_mutes;
DROP TRIGGER set_expired_active_voice_mutes ON active_voice_mutes;
DROP TABLE temporary_rooms;
ALTER TABLE vegans ADD COLUMN notes TEXT;
ALTER TABLE active_members DROP COLUMN guild_snowflake;
ALTER TABLE video_rooms RENAME TO active_video_only_channels;
ALTER TABLE active_voice_mutes
ALTER COLUMN channel_snowflake DROP NOT NULL;
ALTER TABLE active_voice_mutes
ALTER COLUMN target SET NOT NULL;
INSERT INTO active_voice_mutes (
    channel_snowflake,
    created_at,
    expires_in,
    guild_snowflake,
    member_snowflake,
    reason,
    target,
    updated_at
)
SELECT
    NULL,
    created_at,
    expires_in,
    guild_snowflake,
    member_snowflake,
    reason,
    'server',
    updated_at
FROM active_server_voice_mutes;
ALTER TABLE active_voice_mutes
ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_voice_mutes
ALTER COLUMN channel_snowflake DROP NOT NULL;
ALTER TABLE active_text_mutes
ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE command_aliases 
ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_bans 
ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_caps 
ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE active_flags 
ALTER COLUMN channel_snowflake DROP DEFAULT;
ALTER TABLE moderation_logs
ALTER COLUMN guild_snowflake DROP NOT NULL;
ALTER TABLE streaming
RENAME COLUMN guild_snowflake to source_guild_snowflake;
ALTER TABLE streaming
ALTER COLUMN source_guild_snowflake DROP NOT NULL;
ALTER TABLE streaming
ADD COLUMN target_guild_snowflake BIGINT NOT NULL DEFAULT 801609515391778826;
UPDATE streaming
SET target_guild_snowflake = 1347284827350630591
WHERE target_channel_snowflake = 1390814952285012133;
ALTER TABLE streaming
ALTER COLUMN target_guild_snowflake DROP DEFAULT;
ALTER TABLE streaming
RENAME COLUMN target_guild_snowflake TO guild_snowflake;
ALTER TABLE streaming
RENAME COLUMN target_channel_snowflake TO channel_snowflake;
CREATE TABLE permission_entries (
    channel_snowflake bigint,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint,
    group_alias TEXT NOT NULL,
    member_snowflake bigint NOT NULL,
    role_snowflakes bigint[] CONSTRAINT administrators_role_snowflake_not_null NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);
ALTER TABLE permission_entries
ADD CONSTRAINT permission_entries_unique
UNIQUE (
    member_snowflake,
    group_alias,
    guild_snowflake,
    channel_snowflake
);
BEGIN;

INSERT INTO permission_entries (
    channel_snowflake,
    created_at,
    guild_snowflake,
    group_alias,
    member_snowflake,
    updated_at,
    role_snowflakes
)
SELECT
    NULL,
    created_at,
    NULL,
    'developer',
    member_snowflake,
    updated_at,
    ARRAY[]::BIGINT[]
FROM developers;

INSERT INTO permission_entries (
    channel_snowflake,
    created_at,
    guild_snowflake,
    group_alias,
    member_snowflake,
    updated_at,
    role_snowflakes
)
SELECT
    NULL,
    created_at,
    guild_snowflake,
    'administrator',
    member_snowflake,
    updated_at,
    role_snowflakes
FROM administrators;

INSERT INTO permission_entries (
    channel_snowflake,
    created_at,
    guild_snowflake,
    group_alias,
    member_snowflake,
    updated_at,
    role_snowflakes
)
SELECT
    channel_snowflake,
    created_at,
    guild_snowflake,
    'coordinator',
    member_snowflake,
    updated_at,
    ARRAY[]::BIGINT[]
FROM coordinators;

INSERT INTO permission_entries (
    channel_snowflake,
    created_at,
    guild_snowflake,
    group_alias,
    member_snowflake,
    updated_at,
    role_snowflakes
)
SELECT
    channel_snowflake,
    created_at,
    guild_snowflake,
    'moderator',
    member_snowflake,
    updated_at,
    ARRAY[]::BIGINT[]
FROM moderators;
COMMIT;
DROP TABLE moderators;
DROP TABLE administrators;
DROP TABLE developers;
DROP TABLE coordinators;
ALTER TABLE active_video_only_channels
ADD COLUMN expires_in TIMESTAMP WITH TIME ZONE NULL;
ALTER TABLE vegans
ALTER COLUMN created_at DROP NOT NULL;
UPDATE vegans
SET notes = 'No notes provided.'
WHERE notes is NULL;
ALTER TABLE vegans
ALTER COLUMN notes SET NOT NULL;
ALTER TABLE command_aliases
ALTER COLUMN channel_snowflake SET NOT NULL;
ALTER TABLE active_text_mutes
ALTER COLUMN reason SET NOT NULL;
UPDATE active_voice_mutes
SET reason = 'No reason provided.'
WHERE reason is NULL;
ALTER TABLE active_voice_mutes
ALTER COLUMN reason SET NOT NULL;
UPDATE active_bans
SET reason = 'No reason provided.'
WHERE reason is NULL;
ALTER TABLE active_bans
ALTER COLUMN reason SET NOT NULL;
UPDATE active_flags
SET reason = 'No reason provided.'
WHERE reason is NULL;
ALTER TABLE active_flags
ALTER COLUMN reason SET NOT NULL;
ALTER TABLE moderation_logs
ADD COLUMN id SERIAL PRIMARY KEY;
ALTER TABLE active_members
ADD PRIMARY KEY (member_snowflake);

ALTER TABLE administrator_roles
ADD COLUMN group_alias TEXT NOT NULL DEFAULT 'administrator';

ALTER TABLE administrator_roles
ALTER COLUMN group_alias DROP DEFAULT;

ALTER TABLE administrator_roles RENAME TO autoassign_roles;
ALTER TABLE autoassign_roles ADD COLUMN channel_snowflake BIGINT;
ALTER TABLE autoassign_roles RENAME CONSTRAINT administrator_roles_pkey TO autoassign_roles_pkey;
ALTER TABLE roles DROP COLUMN member_snowflake;
INSERT INTO roles (guild_snowflake, channel_snowflake, role_snowflake)
SELECT guild_snowflake, channel_snowflake, role_snowflake
FROM command_aliases
WHERE role_snowflake IS NOT NULL;
ALTER TABLE active_automute_channels RENAME CONSTRAINT active_stages_pkey TO active_automute_channels_pkey;
DROP TABLE active_server_voice_mutes;
ALTER TABLE active_video_only_channels RENAME CONSTRAINT video_rooms_pkey TO active_video_only_channels_pkey;
ALTER SEQUENCE history_id_seq RENAME TO streaming_id_seq;
DELETE FROM moderation_logs WHERE reason='Right-click voice unmute';
