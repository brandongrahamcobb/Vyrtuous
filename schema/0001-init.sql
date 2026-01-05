CREATE USER vyrtuous WITH PASSWORD 'password';
CREATE ROLE
CREATE DATABASE vyrtuous;
CREATE DATABASE
ALTER DATABASE vyrtuous OWNER to vyrtuous;

CREATE TABLE users (
    discord_snowflake BIGINT PRIMARY KEY,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE active_server_voice_mutes (
    guild_id BIGINT NOT NULL,
    discord_snowflake BIGINT NOT NULL,
    expires_in TIMESTAMPTZ,
    reason TEXT,
    PRIMARY KEY (guild_id, discord_snowflake)
);

CREATE TABLE moderation_logs (
    action_type VARCHAR(255) NOT NULL,
    channel_members_voice_count INTEGER DEFAULT 0 NOT NULL,
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL,
    executor_member_snowflake BIGINT NOT NULL,
    expires_at TIMESTAMP,
    guild_members_offline_and_online_member_count INTEGER DEFAULT 0 NOT NULL,
    guild_members_online_count INTEGER DEFAULT 0 NOT NULL,
    guild_members_voice_count INTEGER DEFAULT 0 NOT NULL,
    guild_snowflake BIGINT NOT NULL,
    highest_role VARCHAR(255),
    is_modification BOOLEAN DEFAULT FALSE NOT NULL,
    target_member_snowflake BIGINT,
    reason TEXT,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL
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
    channel_snowflake BIGINT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    owner_snowflake BIGINT NOT NULL,
    room_name TEXT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, room_name)
);

CREATE TABLE moderators (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

CREATE TABLE coordinators (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

CREATE TABLE developers (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, member_snowflake)
);

CREATE TABLE administrators (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    role_snowflakes BIGINT[] NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, member_snowflake)
);

CREATE TABLE active_server_voice_mutes (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    expires_in TIMESTAMPTZ,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, member_snowflake)
);

CREATE TABLE active_text_mutes (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_in TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

CREATE TABLE active_bans (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_in TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

CREATE TABLE active_voice_mutes (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_in TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    target TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

CREATE TABLE temporary_rooms (
    channel_snowflake BIGINT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    room_name TEXT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, room_name)
);

CREATE TABLE active_stages (
    channel_snowflake BIGINT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_in TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake)
);

CREATE TABLE history (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    enabled BOOLEAN DEFAULT FALSE,
    entry_type TEXT NOT NULL,
    guild_snowflake BIGINT NOT NULL,
    id BIGSERIAL PRIMARY KEY,
    snowflakes BIGINT[],
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE vegans (
    channel_snowflake BIGINT NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

CREATE TABLE active_flags (
    channel_snowflake BIGINT NOT NULL DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    expires_in TIMESTAMPTZ,
    guild_snowflake BIGINT NOT NULL,
    member_snowflake BIGINT NOT NULL,
    reason TEXT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake)
);

CREATE TABLE command_aliases (
    alias_type        TEXT NOT NULL CHECK (alias_type IN (
        'vegan', 'carnist', 'voice_mute', 'unvoice_mute', 'ban', 'unban', 'flag', 'unflag', 'text_mute', 'untext_mute', 'role', 'unrole'
    )),
    alias_name        TEXT NOT NULL,
    channel_snowflake BIGINT DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake   BIGINT NOT NULL,
    role_snowflake    BIGINT,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (alias_name, alias_type, guild_snowflake)
);

CREATE TABLE active_caps (
    channel_snowflake BIGINT DEFAULT -1,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    duration_seconds INTEGER NOT NULL,
    guild_snowflake BIGINT NOT NULL,
    moderation_type TEXT NOT NULL CHECK (moderation_type IN ('ban', 'voice_mute', 'text_mute')),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake, moderation_type)
);

CREATE TABLE video_rooms (
    channel_snowflake BIGINT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (channel_snowflake, guild_snowflake)
);

CREATE TABLE administrator_roles (
    created_at TIMESTAMPTZ DEFAULT NOW(),
    guild_snowflake BIGINT NOT NULL,
    role_snowflake BIGINT NOT NULL,
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (guild_snowflake, role_snowflake)
);

CREATE TABLE developer_logs (
    channel_snowflake BIGINT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    developer_snowflakes BIGINT[],
    guild_snowflake BIGINT NOT NULL,
    id UUID PRIMARY KEY,
    message_snowflake BIGINT NOT NULL,
    notes TEXT,
    resolved BOOLEAN NOT NULL DEFAULT FALSE,
    updated_at TIMESTAMPTZ DEFAULT NOW()
);
