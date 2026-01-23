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
DROP TABLE temporary_room_aliases;
DROP TABLE statistic_channels_old;
DROP TABLE stage_coordinators;
DROP TABLE role_permissions;
DROP TABLE moderation_logs_old;
DROP TABLE log_channels;
DROP TABLE active_caps_old;
ALTER TABLE history RENAME TO streaming;
BEGIN;

UPDATE command_aliases
SET alias_type = CASE alias_type
    WHEN 'voice_mute' THEN 'vmute'
    WHEN 'unvoice_mute' THEN 'unvmute'
    WHEN 'text_mute' THEN 'tmute'
    WHEN 'untext_mute' THEN 'untmute'
    ELSE alias_type
END;

COMMIT;

BEGIN;

ALTER TABLE command_aliases
DROP CONSTRAINT command_aliases_alias_type_check;

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
ALTER TABLE your_table_name
RENAME COLUMN developer_snowflakes TO member_snowflakes;
ALTER TABLE active_caps
DROP CONSTRAINT active_caps_moderation_type_check;

UPDATE active_caps
SET moderation_type='vmute'
WHERE moderation_type='voice_mute';

UPDATE active_caps
SET moderation_type='tmute'
WHERE moderation_type='text_mute';
ALTER TABLE active_caps
ADD CONSTRAINT active_caps_moderation_type_check
CHECK (moderation_type IN ('ban', 'vmute', 'tmute'));
ALTER TABLE developer_logs RENAME TO bug_tracking;
CREATE TABLE sysadmins (
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
