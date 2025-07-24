ALTER DATABASE vyrtuous OWNER TO spawd;
DROP TABLE IF EXISTS users;
DROP TABLE IF EXISTS command_aliases;
DROP TABLE IF EXISTS mute_reasons;
DROP TABLE IF EXISTS active_mutes;

CREATE TABLE IF NOT EXISTS users (
    user_id BIGINT PRIMARY KEY,
    mute_channel_ids BIGINT[],              -- an array of channel IDs where user is muted
    manual_mute_channels BIGINT[],
    role_ids BIGINT[],                      -- an array of role IDs associated with the user
    moderator_ids BIGINT[],                 -- an array of channel IDs where user is a moderator
    developer_guild_ids BIGINT[],           -- an array of guild IDs where user has developer rights
    flagged BOOLEAN DEFAULT FALSE,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);
CREATE TABLE command_aliases (
    guild_id BIGINT NOT NULL,
    alias_type TEXT NOT NULL CHECK (alias_type IN ('mute', 'unmute')),
    alias_name TEXT NOT NULL,
    channel_id BIGINT NOT NULL,
    PRIMARY KEY (guild_id, alias_type, alias_name)
);
CREATE TABLE mute_reasons (
    guild_id   BIGINT NOT NULL,
    user_id    BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    reason     TEXT,
    PRIMARY KEY (guild_id, user_id, channel_id)
);
CREATE TABLE active_mutes (
    user_id BIGINT,
    channel_id BIGINT,
    source TEXT CHECK (source IN ('bot', 'manual')),
    PRIMARY KEY (user_id, channel_id)
)

GRANT ALL PRIVILEGES ON DATABASE vyrtuous TO spawd;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO spawd;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO spawd;
GRANT ALL PRIVILEGES ON ALL FUNCTIONS IN SCHEMA public TO spawd;
GRANT USAGE ON SCHEMA public TO spawd;
ALTER DEFAULT PRIVILEGES IN SCHEMA public
GRANT ALL ON TABLES TO spawd;

ALTER DEFAULT PRIVILEGES IN SCHEMA public
GRANT ALL ON SEQUENCES TO spawd;

ALTER DEFAULT PRIVILEGES IN SCHEMA public
GRANT ALL ON FUNCTIONS TO spawd;




