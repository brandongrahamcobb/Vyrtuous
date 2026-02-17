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

ALTER TABLE your_table
ADD CONSTRAINT unique_target_source UNIQUE (target_channel_snowflake, source_channel_snowflake);
