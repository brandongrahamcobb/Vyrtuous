CREATE DATABASE vyrtuous;
ALTER DATABASE vyrtuous OWNER TO spawd;
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
    user_id BIGINT PRIMARY KEY,
    ban_channel_ids BIGINT[],
    mute_channel_ids BIGINT[],
    manual_mute_channels BIGINT[],
    moderator_ids BIGINT[],
    moderator_channel_ids BIGINT[],
    coordinator_ids BIGINT[],
    coordinator_channel_ids BIGINT[],
    developer_guild_ids BIGINT[],
    flagged_channel_ids BIGINT[],
    going_vegan_channel_ids BIGINT[],
    server_mute_guild_ids BIGINT[],
    server_muter_guild_ids BIGINT[],
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    created_at TIMESTAMPTZ DEFAULT NOW()
);
CREATE TABLE text_mutes (
    guild_id BIGINT NOT NULL,
    user_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    issuer_id BIGINT NOT NULL,
    reason TEXT NOT NULL,
    source TEXT NOT NULL,
    expires_at TIMESTAMPTZ,
    PRIMARY KEY (guild_id, user_id, channel_id)
);
-- CREATE TABLE text_mutes (
   -- user_id BIGINT NOT NULL,
   -- channel_id BIGINT NOT NULL,
   -- issuer_id BIGINT NOT NULL,
   -- reason TEXT NOT NULL,
   -- source TEXT NOT NULL,
   -- expires_at TIMESTAMPTZ,
   -- PRIMARY KEY (user_id, channel_id)
--);

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

CREATE TABLE IF NOT EXISTS mute_reasons (
    guild_id BIGINT NOT NULL,
    user_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    reason TEXT,
    PRIMARY KEY (guild_id, user_id, channel_id)
);


CREATE TABLE IF NOT EXISTS ban_reasons (
    guild_id BIGINT NOT NULL,
    user_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    reason TEXT,
    PRIMARY KEY (guild_id, user_id, channel_id)
);

-- Currently active bans
CREATE TABLE IF NOT EXISTS active_bans (
    user_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    PRIMARY KEY (user_id, channel_id)
);

-- Currently active mutes

CREATE TABLE IF NOT EXISTS active_mutes (
    user_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    source TEXT CHECK (source IN ('bot', 'bot_owner', 'manual', 'owner')),
    issuer_id BIGINT,
    expires_at TIMESTAMPTZ,
    PRIMARY KEY (user_id, channel_id)
);


-- Expiring bans
CREATE TABLE IF NOT EXISTS ban_expirations (
    user_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ NOT NULL,
    PRIMARY KEY (user_id, channel_id)
);

CREATE TABLE IF NOT EXISTS server_mute_reasons (
    guild_id BIGINT NOT NULL,
    user_id BIGINT NOT NULL,
    reason TEXT,
    PRIMARY KEY (guild_id, user_id)
);
CREATE TABLE moderation_logs (
    id BIGSERIAL PRIMARY KEY,
    action_type TEXT NOT NULL,
    target_user_id BIGINT,
    executor_user_id BIGINT NOT NULL,
    guild_id BIGINT NOT NULL,
    channel_id BIGINT,
    reason TEXT,
    metadata JSONB DEFAULT '{}'::jsonb,                
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

CREATE TABLE mute_reasons (
    guild_id BIGINT NOT NULL,
    user_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    reason TEXT NOT NULL,
    PRIMARY KEY (guild_id, user_id, channel_id)
);

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




