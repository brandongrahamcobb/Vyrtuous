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
ALTER TABLE administrators DROP COLUMN display_name;
DROP TABLE ban_roles;
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
DROP TABLE roles;
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
ALTER TABLE ONLY active_voice_mutes
DROP CONSTRAINT active_voice_mutes_pkey;
ALTER TABLE ONLY active_voice_mutes
ADD CONSTRAINT active_voice_mutes_pkey
PRIMARY KEY (guild_snowflake, member_snowflake, target);
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
ADD COLUMN target_guild_snowflake BIGINT NOT NULL;
ALTER TABLE streaming
RENAME COLUMN target_guild_snowflake TO guild_snowflake;
ALTER TABLE streaming
RENAME COLUMN target_channel_snowflake TO channel_snowflake;

CREATE TABLE permission_levels (
    channel_snowflake bigint,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint,
    level_name TEXT NOT NULL,
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
);
