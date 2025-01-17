-- schema.sql
-- This script drops existing tables and creates a new schema for the TagManager.

-- Drop existing tables if they exist
DROP TABLE IF EXISTS loop_configs;
DROP TABLE IF EXISTS tags;

-- Create the 'tags' table with the updated schema
CREATE TABLE IF NOT EXISTS tags (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    location_id BIGINT NOT NULL, -- Typically the Guild ID
    owner_id BIGINT NOT NULL,     -- Discord User ID
    content TEXT,
    attachment_url TEXT,
    tag_type VARCHAR(50) DEFAULT 'default',
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(name, location_id, owner_id)
);

-- Create the 'loop_configs' table
CREATE TABLE IF NOT EXISTS loop_configs (
    guild_id BIGINT PRIMARY KEY,
    channel_id BIGINT,
    enabled BOOLEAN NOT NULL DEFAULT TRUE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Optional: Create indexes to improve query performance
CREATE INDEX IF NOT EXISTS idx_tags_name ON tags(name);
CREATE INDEX IF NOT EXISTS idx_tags_location ON tags(location_id);
CREATE INDEX IF NOT EXISTS idx_loop_configs_enabled ON loop_configs(enabled);
