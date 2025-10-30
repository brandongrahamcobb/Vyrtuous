ALTER DATABASE vyrtuous OWNER TO vyrtuous;
ALTER ROLE vyrtuous WITH SUPERUSER;
DROP TABLE IF EXISTS users;
DROP TABLE IF EXISTS command_aliases;
DROP TABLE IF EXISTS mute_reasons;
DROP TABLE IF EXISTS active_mutes;
DROP TABLE IF EXISTS ban_reasons;
DROP TABLE IF EXISTS active_bans;
DROP TABLE IF EXISTS ban_roles;
DROP TABLE IF EXISTS mute_roles;
DROP TABLE IF EXISTS moderation_logs;
DROP TABLE IF EXISTS mute_reasons;
DROP TABLE IF EXISTS server_mute_reasons;
DROP TABLE IF EXISTS server_muter_ids;
DROP TABLE IF EXISTS ban_expirations;
DROP TABLE IF EXISTS text_mutes;

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
CREATE TABLE active_text_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    reason TEXT,
    expires_at TIMESTAMPTZ,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

CREATE TABLE active_caps (
    guild_id   BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    moderation_type   TEXT   NOT NULL,  -- 'ban', 'mute', or 'tmute'
    duration   TEXT   NOT NULL,  -- parsed duration string, e.g. '24h'
    PRIMARY KEY (guild_id, channel_id, moderation_type)
);
CREATE TABLE IF NOT EXISTS command_aliases (
    guild_id BIGINT NOT NULL,
    alias_type TEXT NOT NULL CHECK (alias_type IN (
        'cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute', 'role', 'unrole'
    )),
    alias_name TEXT NOT NULL,
    channel_id BIGINT,
    role_id BIGINT,
    PRIMARY KEY (guild_id, alias_type, alias_name)
);

-- Currently active bans
CREATE TABLE IF NOT EXISTS active_bans (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

CREATE TABLE IF NOT EXISTS active_server_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake)
);

CREATE TABLE IF NOT EXISTS active_cows (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

CREATE TABLE IF NOT EXISTS active_flags (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

CREATE TABLE IF NOT EXISTS active_stages (
    guild_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    initiator_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    PRIMARY KEY (guild_id, channel_id)
);

CREATE TABLE IF NOT EXISTS stage_coordinators (
    guild_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    PRIMARY KEY (guild_id, channel_id, discord_snowflake)
);

CREATE TABLE IF NOT EXISTS active_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    target TEXT CHECK (target IN ('room', 'user')) DEFAULT 'user',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, target)
);

CREATE TABLE moderation_logs (
    id BIGSERIAL PRIMARY KEY,
    action_type TEXT NOT NULL,
    target_discord_snowflake BIGINT,
    executor_discord_snowflake BIGINT NOT NULL,
    guild_id BIGINT NOT NULL,
    channel_id BIGINT,
    reason TEXT,
    metadata JSONB DEFAULT '{}'::jsonb,                
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);
CREATE TABLE IF NOT EXISTS log_channels (
    guild_id   BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    type       TEXT DEFAULT 'general',
    snowflakes BIGINT[],
    enabled BOOLEAN DEFAULT FALSE,
    PRIMARY KEY (guild_id, channel_id)
);

GRANT ALL PRIVILEGES ON DATABASE vyrtuous TO vyrtuous;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO vyrtuous;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO vyrtuous;
GRANT ALL PRIVILEGES ON ALL FUNCTIONS IN SCHEMA public TO vyrtuous;
GRANT USAGE ON SCHEMA public TO vyrtuous;
ALTER DEFAULT PRIVILEGES IN SCHEMA public
GRANT ALL ON TABLES TO vyrtuous;

ALTER DEFAULT PRIVILEGES IN SCHEMA public
GRANT ALL ON SEQUENCES TO vyrtuous;

ALTER DEFAULT PRIVILEGES IN SCHEMA public
GRANT ALL ON FUNCTIONS TO vyrtuous;




