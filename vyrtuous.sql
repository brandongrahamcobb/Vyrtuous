DROP TABLE IF EXISTS users;
DROP TABLE IF EXISTS command_aliases;
CREATE TABLE IF NOT EXISTS users (
    user_id BIGINT PRIMARY KEY,
    mute_channel_ids BIGINT[],              -- an array of channel IDs where user is muted
    role_ids BIGINT[],                      -- an array of role IDs associated with the user
    moderator_ids BIGINT[],                 -- an array of channel IDs where user is a moderator
    developer_guild_ids BIGINT[],           -- an array of guild IDs where user has developer rights
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
