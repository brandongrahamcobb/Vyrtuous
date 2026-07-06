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
ALTER TABLE active_voice_mutes ALTER COLUMN channel_snowflake SET DEFAULT (-1)::bigint;
