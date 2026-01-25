ALTER TABLE active_bans ADD COLUMN expired BOOLEAN DEFAULT FALSE;
ALTER TABLE active_text_mutes ADD COLUMN expired BOOLEAN DEFAULT FALSE;
ALTER TABLE active_voice_mutes ADD COLUMN expired BOOLEAN DEFAULT FALSE;
ALTER TABLE active_stages ADD COLUMN expired BOOLEAN DEFAULT FALSE;
ALTER TABLE developer_logs ADD COLUMN expired BOOLEAN DEFAULT FALSE;
ALTER TABLE developers DROP CONSTRAINT developers_pkey;
BEGIN;

WITH ranked AS (
  SELECT
    *,
    ROW_NUMBER() OVER (PARTITION BY member_snowflake ORDER BY updated_at DESC) AS rn
  FROM developers
)
DELETE FROM developers
WHERE (member_snowflake, updated_at) IN (
  SELECT member_snowflake, updated_at
  FROM ranked
  WHERE rn > 1
);

ALTER TABLE developers DROP COLUMN guild_snowflake;
ALTER TABLE developers ADD PRIMARY KEY (member_snowflake);

COMMIT;
CREATE OR REPLACE FUNCTION set_expired()
RETURNS TRIGGER AS $$
BEGIN
    NEW.expired := (NEW.expires_in IS NOT NULL AND NEW.expires_in < NOW());
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;
CREATE TRIGGER set_expired_active_bans
BEFORE INSERT OR UPDATE ON active_bans
FOR EACH ROW
EXECUTE FUNCTION set_expired();

CREATE TRIGGER set_expired_active_text_mutes
BEFORE INSERT OR UPDATE ON active_text_mutes
FOR EACH ROW
EXECUTE FUNCTION set_expired();

CREATE TRIGGER set_expired_active_voice_mutes
BEFORE INSERT OR UPDATE ON active_voice_mutes
FOR EACH ROW
EXECUTE FUNCTION set_expired();

CREATE TRIGGER set_expired_active_stages
BEFORE INSERT OR UPDATE ON active_stages
FOR EACH ROW
EXECUTE FUNCTION set_expired();

CREATE TRIGGER set_expired_developer_logs
BEFORE INSERT OR UPDATE ON developer_logs
FOR EACH ROW
EXECUTE FUNCTION set_expired();

DROP TABLE temporary_room_aliases;
DROP TABLE statistic_channels_old;
DROP TABLE stage_coordinators;
DROP TABLE role_permissions;
DROP TABLE moderation_logs_old;
DROP TABLE log_channels;
DROP TABLE active_caps_old;
ALTER TABLE history RENAME TO streaming;
BEGIN;

ALTER TABLE command_aliases
DROP CONSTRAINT IF EXISTS command_aliases_alias_type_check;

DELETE FROM command_aliases
WHERE alias_type IS NULL;

UPDATE command_aliases
SET alias_type = CASE alias_type
    WHEN 'voice_mute' THEN 'vmute'
    WHEN 'unvoice_mute' THEN 'unvmute'
    WHEN 'text_mute' THEN 'tmute'
    WHEN 'untext_mute' THEN 'untmute'
    ELSE alias_type
END;

DELETE FROM command_aliases
WHERE alias_type NOT IN (
    'vegan','carnist','vmute','unvmute','ban','unban','flag','unflag','tmute','untmute','role','unrole'
);

ALTER TABLE command_aliases
ADD CONSTRAINT command_aliases_alias_type_check
CHECK (alias_type IN (
    'vegan','carnist','vmute','unvmute','ban','unban','flag','unflag','tmute','untmute','role','unrole'
));

COMMIT;

ALTER TABLE vegans
DROP CONSTRAINT vegans_pkey,
DROP COLUMN channel_snowflake,
ADD PRIMARY KEY (guild_snowflake, member_snowflake);
ALTER TABLE active_stages DROP COLUMN member_snowflake;
ALTER TABLE active_bans
ADD COLUMN role_snowflake BIGINT;
ALTER TABLE active_text_mutes
ADD COLUMN role_snowflake BIGINT;
DELETE FROM command_aliases
WHERE alias_type IN ('unban', 'unrole', 'untmute', 'unvmute', 'unflag', 'carnist');
ALTER TABLE developer_logs
RENAME COLUMN developer_snowflakes TO member_snowflakes;
ALTER TABLE developer_logs RENAME TO bug_tracking;
CREATE TABLE sysadmin (
    id BOOLEAN PRIMARY KEY DEFAULT TRUE,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    CHECK (id)
);
CREATE TABLE guild_owners (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, member_snowflake)
);
ALTER TABLE active_bans RENAME TO active_hides;
ALTER TABLE active_hides RENAME CONSTRAINT active_bans_pkey1 TO active_hides_pkey;
CREATE TABLE active_bans (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expired BOOLEAN DEFAULT FALSE,
    expires_in TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    role_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);
CREATE OR REPLACE FUNCTION set_expired()
RETURNS TRIGGER AS $$
BEGIN
    NEW.expired := (NEW.expires_in IS NOT NULL AND NEW.expires_in < NOW());
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;
CREATE TRIGGER set_expired_active_bans
BEFORE INSERT OR UPDATE ON active_bans
FOR EACH ROW
EXECUTE FUNCTION set_expired();
ALTER TABLE bug_tracking RENAME CONSTRAINT developer_logs_pkey TO bug_tracking_pkey;
ALTER TABLE active_voice_mutes RENAME CONSTRAINT active_voice_mutes_pkey1 TO active_voice_mutes_pkey;
ALTER TABLE streaming RENAME CONSTRAINT history_pkey TO streaming_pkey;
ALTER TABLE temporary_rooms RENAME CONSTRAINT temporary_rooms_pkey1 TO temporary_rooms_pkey;

