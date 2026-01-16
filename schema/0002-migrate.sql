ALTER TABLE active_bans
ADD COLUMN expired boolean GENERATED ALWAYS AS (
    expires_in IS NOT NULL AND expires_in < NOW()
) STORED;
ALTER TABLE active_text_mutes
ADD COLUMN expired boolean GENERATED ALWAYS AS (
    expires_in IS NOT NULL AND expires_in < NOW()
) STORED;
ALTER TABLE active_voice_mutes
ADD COLUMN expired boolean GENERATED ALWAYS AS (
    expires_in IS NOT NULL AND expires_in < NOW()
) STORED;
ALTER TABLE active_stages
ADD COLUMN expired boolean GENERATED ALWAYS AS (
    expires_in IS NOT NULL AND expires_in < NOW()
) STORED;
ALTER TABLE developer_logs
ADD COLUMN expired boolean GENERATED ALWAYS AS (
    expires_in IS NOT NULL AND expires_in < NOW()
) STORED;
ALTER TABLE developers DROP CONSTRAINT developers_pkey;
ALTER TABLE developers DROP COLUMN guild_snowflake;
ALTER TABLE developers ADD PRIMARY KEY (member_snowflake);
