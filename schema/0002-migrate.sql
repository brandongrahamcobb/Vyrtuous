ALTER TABLE public.streaming
DROP COLUMN IF EXISTS channel_snowflake,
DROP COLUMN IF EXISTS enabled,
DROP COLUMN IF EXISTS entry_type,
DROP COLUMN IF EXISTS snowflakes;

ALTER TABLE public.streaming
ADD COLUMN IF NOT EXISTS target_channel_snowflake bigint NOT NULL,
ADD COLUMN IF NOT EXISTS source_channel_snowflake bigint NOT NULL;

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
