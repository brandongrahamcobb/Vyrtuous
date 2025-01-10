-- SQL File to Set Up Tag System Database

-- Drop existing tables if they exist
-- DROP TABLE IF EXISTS tags;
-- DROP TABLE IF EXISTS loop_configs;

-- Create the tags table
-- CREATE TABLE tags (
--    id SERIAL PRIMARY KEY,
--    name VARCHAR(100) NOT NULL,
--    location_id BIGINT NOT NULL,
--    content TEXT NULL,
--    attachment_url TEXT NULL,
--    owner_id BIGINT NOT NULL,
--    tag_type VARCHAR(10) NOT NULL CHECK (tag_type IN ('default', 'loop')),
--    UNIQUE (name, location_id, owner_id)
-- );

CREATE TABLE tags (
    id SERIAL PRIMARY KEY,
    name VARCHAR(100) NOT NULL,
    location_id BIGINT NOT NULL, -- Guild ID
    content TEXT,
    attachment_url TEXT,
    owner_id BIGINT NOT NULL, -- User ID
    tag_type VARCHAR(50) DEFAULT 'default',
    UNIQUE (LOWER(name), location_id, owner_id)
);

-- Create the loop_configs table
CREATE TABLE loop_configs (
    guild_id BIGINT PRIMARY KEY,
    channel_id BIGINT NOT NULL,
    enabled BOOLEAN NOT NULL DEFAULT FALSE
);
