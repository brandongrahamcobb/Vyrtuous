ALTER TABLE public.streaming
DROP COLUMN IF EXISTS enabled,
DROP COLUMN IF EXISTS entry_type,
DROP COLUMN IF EXISTS snowflakes;

ALTER TABLE public.streaming
RENAME COLUMN channel_snowflake TO target_channel_snowflake;
ALTER TABLE public.streaming
ALTER COLUMN target_channel_snowflake SET NOT NULL;
ALTER TABLE public.streaming
ADD COLUMN IF NOT EXISTS source_channel_snowflake bigint;

ALTER TABLE public.streaming
DROP CONSTRAINT IF EXISTS streaming_pkey;

ALTER TABLE streaming
ADD CONSTRAINT unique_target_source UNIQUE (target_channel_snowflake, source_channel_snowflake);

ALTER TABLE moderation_logs RENAME COLUMN channel_members_voice_count TO current_channel_members;
ALTER TABLE moderation_logs RENAME COLUMN guild_members_offline_and_online_member_count TO total_guild_members;
ALTER TABLE moderation_logs RENAME COLUMN guild_members_online_count TO online_members;
ALTER TABLE moderation_logs RENAME COLUMN guild_members_voice_count TO total_voice_members;
ALTER TABLE moderation_logs RENAME COLUMN infraction_type TO identifier;
ALTER TABLE moderation_logs RENAME COLUMN executor_member_snowflake TO author_snowflake;
ALTER TABLE moderation_logs RENAME COLUMN target_member_snowflake TO target_snowflake;
CREATE TABLE uploads (
    command_name TEXT NOT NULL,
    file_bytes BYTEA NOT NULL,
    filename TEXT NOT NULL,
    tag TEXT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (command_name, tag)
);
ALTER TABLE streaming
ALTER COLUMN source_channel_snowflake DROP NOT NULL;
ALTER TABLE active_bans
ADD COLUMN display_name TEXT;
ALTER TABLE active_bans
ADD COLUMN blacklisted BOOLEAN;
ALTER TABLE guild_owners
ADD COLUMN display_name TEXT;
ALTER TABLE active_voice_mutes
ADD COLUMN display_name TEXT;
ALTER TABLE active_text_mutes
ADD COLUMN display_name TEXT;
ALTER TABLE developers
ADD COLUMN display_name TEXT;
ALTER TABLE moderators
ADD COLUMN display_name TEXT;
ALTER TABLE coordinators
ADD COLUMN display_name TEXT;
ALTER TABLE sysadmin
ADD COLUMN display_name TEXT;
ALTER TABLE administrators
ADD COLUMN display_name TEXT;
ALTER TABLE active_flags
ADD COLUMN display_name TEXT;
ALTER TABLE vegans
ADD COLUMN display_name TEXT;
ALTER TABLE active_server_voice_mutes
ADD COLUMN display_name TEXT;


CREATE TABLE active_members (
    created_at timestamp with time zone DEFAULT now(),
    display_name TEXT,
    guild_snowflake bigint NOT NULL,
    last_active timestamp with time zone DEFAULT now(),
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);

