-- 1. Create the new table if it doesn't exist
DROP TABLE ban_expirations;
DROP TABLE channel_roles;
CREATE TABLE IF NOT EXISTS active_bans_new (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

-- 2. Copy data from old table to new table
INSERT INTO active_bans_new (guild_id, discord_snowflake, channel_id, expires_at, reason)
SELECT 801609515391778826, user_id, channel_id, expires_at, NULL
FROM active_bans;

-- 3. (Optional) Rename old table for backup
ALTER TABLE active_bans RENAME TO active_bans_old;

-- 4. Rename new table to the original name
ALTER TABLE active_bans_new RENAME TO active_bans;

-- 5. (Optional) Check your data
SELECT * FROM active_bans LIMIT 10;
DROP TABLE active_bans_old;

-- 1. Create the new table if it doesn't exist
CREATE TABLE IF NOT EXISTS active_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

-- 2. Copy data from old table, merging mute_reasons
INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason)
SELECT
    801609515391778826 AS guild_id,
    am.user_id AS discord_snowflake,
    am.channel_id,
    am.expires_at,
    mr.reason
FROM active_mutes am
LEFT JOIN mute_reasons mr
    ON mr.guild_id = 801609515391778826
    AND mr.user_id = am.user_id
    AND mr.channel_id = am.channel_id;

-- 3. (Optional) Backup old table
ALTER TABLE active_mutes RENAME TO active_mutes_old;

-- 4. (Optional) Verify migrated data
SELECT * FROM active_voice_mutes LIMIT 10;
DROP TABLE active_mutes_old;
DROP TABLE mute_reasons;
-- 1. Create the new table if it doesn't exist
CREATE TABLE IF NOT EXISTS active_server_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake)
);

-- 2. Copy data from old table
-- No expires_at exists in server_mute_reasons, so it will be NULL
INSERT INTO active_server_voice_mutes (guild_id, discord_snowflake, expires_at, reason)
SELECT
    guild_id,
    user_id AS discord_snowflake,
    NULL AS expires_at,
    reason
FROM server_mute_reasons;

-- 3. (Optional) Backup old table
ALTER TABLE server_mute_reasons RENAME TO server_mute_reasons_old;

-- 4. (Optional) Verify migrated data
SELECT * FROM active_server_voice_mutes LIMIT 10;
DROP TABLE server_mute_reasons_old;
-- 1. Create the new table if it doesn't exist
CREATE TABLE IF NOT EXISTS active_text_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    reason TEXT,
    expires_at TIMESTAMPTZ,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

-- 2. Copy data from old table
INSERT INTO active_text_mutes (guild_id, discord_snowflake, channel_id, reason, expires_at)
SELECT
    guild_id,
    user_id AS discord_snowflake,
    channel_id,
    reason,
    expires_at
FROM text_mutes;

-- 3. (Optional) Backup old table
ALTER TABLE text_mutes RENAME TO text_mutes_old;
DROP TABLE text_mutes_old;

-- 4. (Optional) Verify migrated data
SELECT * FROM active_text_mutes LIMIT 10;
UPDATE active_bans ab
SET reason = br.reason
FROM ban_reasons br
WHERE ab.guild_id = br.guild_id
  AND ab.discord_snowflake = br.user_id
  AND ab.channel_id = br.channel_id
  AND ab.reason IS NULL;

-- 2. Insert any ban_reasons rows that don’t already exist in active_bans
INSERT INTO active_bans (guild_id, discord_snowflake, channel_id, reason)
SELECT br.guild_id, br.user_id, br.channel_id, br.reason
FROM ban_reasons br
LEFT JOIN active_bans ab
  ON ab.guild_id = br.guild_id
 AND ab.discord_snowflake = br.user_id
 AND ab.channel_id = br.channel_id
WHERE ab.guild_id IS NULL;

DROP TABLE ban_reasons;

ALTER TABLE moderation_logs RENAME TO moderation_logs_old;

-- 2. Create new schema (use BIGINT instead of BIGSERIAL so we control the sequence)
CREATE TABLE moderation_logs (
    id BIGINT PRIMARY KEY,
    action_type TEXT NOT NULL,
    target_discord_snowflake BIGINT,
    executor_discord_snowflake BIGINT NOT NULL,
    guild_id BIGINT NOT NULL,
    channel_id BIGINT,
    reason TEXT,
    metadata JSONB DEFAULT '{}'::jsonb,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- 3. Copy data over
INSERT INTO moderation_logs (
    id,
    action_type,
    target_discord_snowflake,
    executor_discord_snowflake,
    guild_id,
    channel_id,
    reason,
    metadata,
    created_at
)
SELECT
    id,
    action_type,
    target_user_id,
    executor_user_id,
    guild_id,
    channel_id,
    reason,
    metadata,
    created_at
FROM moderation_logs_old;

ALTER TABLE moderation_logs ALTER COLUMN id DROP DEFAULT;

-- 5. Drop the leftover auto-created sequence if it exists
DROP SEQUENCE IF EXISTS moderation_logs_id_seq1;

-- 6. Recreate the original sequence
CREATE SEQUENCE moderation_logs_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;

-- 7. Attach the sequence to the id column
ALTER TABLE moderation_logs ALTER COLUMN id SET DEFAULT nextval('moderation_logs_id_seq');

-- 8. Reset the sequence to the current max id
SELECT setval('moderation_logs_id_seq', (SELECT COALESCE(MAX(id), 1) FROM moderation_logs), true);

-- 9. Make the sequence owned by the column
ALTER SEQUENCE moderation_logs_id_seq OWNED BY moderation_logs.id;

-- 10. Drop the old backup table now that everything is safe
DROP TABLE moderation_logs_old;

CREATE TABLE IF NOT EXISTS active_flags (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

INSERT INTO active_flags (guild_id, discord_snowflake, channel_id, expires_at, reason)
SELECT
    801609515391778826 AS guild_id,
    user_id AS discord_snowflake,
    channel_id,
    NULL AS expires_at,
    NULL AS reason
FROM users,
UNNEST(flagged_channel_ids) AS channel_id
WHERE flagged_channel_ids IS NOT NULL;

CREATE TABLE IF NOT EXISTS active_cows (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

INSERT INTO active_cows (guild_id, discord_snowflake, channel_id)
SELECT
    801609515391778826,
    user_id,
    unnest(array_remove(going_vegan_channel_ids, NULL))
FROM users
WHERE going_vegan_channel_ids IS NOT NULL
  AND array_length(array_remove(going_vegan_channel_ids, NULL), 1) > 0;

ALTER TABLE users RENAME TO users_old;

-- 2. Create new table
CREATE TABLE IF NOT EXISTS users (
    discord_snowflake BIGINT PRIMARY KEY,
    moderator_channel_ids BIGINT[],
    coordinator_channel_ids BIGINT[],
    developer_guild_ids BIGINT[],
    server_mute_guild_ids BIGINT[],
    server_muter_guild_ids BIGINT[],
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- 3. Copy data from old table to new table (only columns that exist in new table)
INSERT INTO users (discord_snowflake, moderator_channel_ids, coordinator_channel_ids, developer_guild_ids, server_mute_guild_ids, server_muter_guild_ids, updated_at, created_at)
SELECT user_id, moderator_channel_ids, coordinator_channel_ids, developer_guild_ids, server_mute_guild_ids, server_muter_guild_ids, updated_at, created_at
FROM users_old;

DROP TABLE users_old;

