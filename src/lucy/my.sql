-- ============================================================
-- 0. Drop Existing Tables (if they exist) for a Fresh Setup
-- ============================================================
DROP TABLE IF EXISTS tags CASCADE;
DROP TABLE IF EXISTS loop_configs CASCADE;

-- ============================================================
-- 1. Create Table: tags
--    This table stores user-defined tags with content and/or attachments.
-- ============================================================
CREATE TABLE IF NOT EXISTS tags (
    id              SERIAL PRIMARY KEY,                -- Unique ID for each tag
    name            TEXT NOT NULL,                     -- Tag name
    location_id     BIGINT NOT NULL,                   -- Guild/Server ID
    content         TEXT,                              -- Tag text content
    attachment_url  TEXT,                              -- Attachment URL (if any)
    owner_id        BIGINT NOT NULL,                   -- User ID of the tag creator
    tag_type        TEXT NOT NULL DEFAULT 'default',   -- Tag type ('default' or 'loop')
    created_at      TIMESTAMP DEFAULT NOW(),           -- Timestamp for creation
    updated_at      TIMESTAMP DEFAULT NOW()            -- Timestamp for last update
);

-- Ensure tag names are unique per guild/location and owner
CREATE UNIQUE INDEX IF NOT EXISTS idx_tags_unique
ON tags (location_id, LOWER(name), owner_id);

-- ============================================================
-- 2. Create Table: loop_configs
--    This table stores configurations for random loop messages.
-- ============================================================
CREATE TABLE IF NOT EXISTS loop_configs (
    guild_id        BIGINT PRIMARY KEY,                -- Guild/Server ID (unique)
    channel_id      BIGINT NOT NULL,                   -- Channel ID where messages are sent
    enabled         BOOLEAN NOT NULL DEFAULT FALSE,    -- Whether the loop is active
    updated_at      TIMESTAMP DEFAULT NOW()            -- Timestamp for the last config update
);

-- ============================================================
-- 3. Add Triggers to Update "updated_at" on Changes
-- ============================================================

-- Update "updated_at" timestamp for tags on update
CREATE OR REPLACE FUNCTION update_tags_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER trigger_update_tags_updated_at
BEFORE UPDATE ON tags
FOR EACH ROW
EXECUTE FUNCTION update_tags_updated_at();

-- Update "updated_at" timestamp for loop_configs on update
CREATE OR REPLACE FUNCTION update_loop_configs_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER trigger_update_loop_configs_updated_at
BEFORE UPDATE ON loop_configs
FOR EACH ROW
EXECUTE FUNCTION update_loop_configs_updated_at();
