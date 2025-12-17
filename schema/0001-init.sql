ALTER ROLE vyrtuous WITH SUPERUSER LOGIN PASSWORD 'password';
ALTER DATABASE vyrtuous OWNER TO vyrtuous;

CREATE TABLE users (
    discord_snowflake BIGINT PRIMARY KEY,
    moderator_channel_ids BIGINT[],
    coordinator_channel_ids BIGINT[],
    developer_guild_ids BIGINT[],
    administrator_guild_ids BIGINT[],
    administrator_role_ids BIGINT[],
    server_mute_guild_ids BIGINT[],
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE active_text_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    reason TEXT,
    expires_at TIMESTAMPTZ,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name)
);

CREATE TABLE active_caps (
    guild_id BIGINT NOT NULL,
    channel_id BIGINT DEFAULT -1,
    room_name TEXT DEFAULT '',
    moderation_type TEXT NOT NULL CHECK (moderation_type IN ('ban', 'mute', 'tmute')),
    duration_seconds INTEGER NOT NULL,
    PRIMARY KEY (guild_id, channel_id, room_name, moderation_type)
);

CREATE TABLE command_aliases (
    guild_id   BIGINT NOT NULL,
    alias_type TEXT NOT NULL CHECK (alias_type IN (
        'cow', 'uncow', 'mute', 'unmute', 'ban', 'unban', 'flag', 'unflag', 'tmute', 'untmute', 'role', 'unrole'
    )),
    alias_name TEXT NOT NULL,
    channel_id BIGINT DEFAULT -1,
    role_id    BIGINT,
    PRIMARY KEY (guild_id, alias_type, alias_name)
);

CREATE TABLE active_bans (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name)
);

CREATE TABLE active_server_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake)
);

CREATE TABLE active_cows (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

CREATE TABLE active_flags (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake, channel_id)
);

CREATE TABLE active_stages (
    guild_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    initiator_id BIGINT NOT NULL,
    expires_at TIMESTAMPTZ,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, channel_id, room_name)
);

CREATE TABLE stage_coordinators (
    guild_id BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    discord_snowflake BIGINT NOT NULL,
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name)
);

CREATE TABLE active_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    channel_id BIGINT NOT NULL DEFAULT -1,
    expires_at TIMESTAMPTZ,
    reason TEXT,
    target TEXT CHECK (target IN ('room', 'user')) DEFAULT 'user',
    room_name TEXT DEFAULT '',
    PRIMARY KEY (guild_id, discord_snowflake, channel_id, room_name, target)
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

CREATE TABLE statistic_channels (
    guild_id   BIGINT NOT NULL,
    channel_id BIGINT NOT NULL,
    type       TEXT DEFAULT 'general',
    snowflakes BIGINT[],
    enabled BOOLEAN DEFAULT FALSE,
    PRIMARY KEY (guild_id, channel_id)
);

CREATE TABLE temporary_rooms (
    guild_snowflake BIGINT NOT NULL,
    room_name TEXT NOT NULL,
    owner_snowflake BIGINT NOT NULL,
    room_snowflake BIGINT,
    PRIMARY KEY (guild_snowflake, room_name)
);
