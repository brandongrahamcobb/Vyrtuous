-- Drop existing tables if needed (Ensure you have backups if necessary)
-- Uncomment the following lines if you need to reset the tables
 DROP TABLE IF EXISTS borrowed_tags CASCADE;
 DROP TABLE IF EXISTS tags CASCADE;
 DROP TABLE IF EXISTS loop_configs CASCADE;
 DROP TABLE IF EXISTS users CASCADE;

-- 1. Create the 'tags' table without the UNIQUE constraint
-- CREATE TABLE IF NOT EXISTS tags (
--     tag_id SERIAL PRIMARY KEY,
--     name VARCHAR(255) NOT NULL,
--     location_id BIGINT NOT NULL, -- Typically the Guild ID
--     owner_id BIGINT NOT NULL,     -- Discord User ID
--     content TEXT,
--     attachment_url TEXT,
--     tag_type VARCHAR(50) DEFAULT 'default',
--     created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
-- );
-- 
-- -- 1.a. Create a unique index for LOWER(name), location_id, owner_id
-- CREATE UNIQUE INDEX IF NOT EXISTS idx_tags_lower_name_location_owner 
-- ON tags (LOWER(name), location_id, owner_id);
-- 
-- -- 2. Create the 'loop_configs' table
-- CREATE TABLE IF NOT EXISTS loop_configs (
--     guild_id BIGINT PRIMARY KEY,
--     channel_id BIGINT,
--     enabled BOOLEAN DEFAULT FALSE,
--     updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
-- );
-- 
-- -- 3. Create the 'borrowed_tags' table
-- CREATE TABLE IF NOT EXISTS borrowed_tags (
--     borrow_id SERIAL PRIMARY KEY,
--     borrower_id BIGINT NOT NULL,       -- Discord User ID who borrowed the tag
--     original_tag_id INT NOT NULL,      -- References the original tag
--     borrowed_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
--     returned_at TIMESTAMP WITH TIME ZONE,
--     FOREIGN KEY (original_tag_id) REFERENCES tags(tag_id) ON DELETE CASCADE
-- );
-- 
-- -- 3.a. Create a partial unique index to ensure one active borrow per user per tag
-- CREATE UNIQUE INDEX IF NOT EXISTS idx_borrowed_tags_unique_active 
-- ON borrowed_tags (borrower_id, original_tag_id) 
-- WHERE returned_at IS NULL;
-- 
-- -- 4. Create the 'users' table
-- CREATE TABLE IF NOT EXISTS users (
--     id BIGINT PRIMARY KEY,
--     name VARCHAR(255) NOT NULL,
--     create_date TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
--     level INT NOT NULL DEFAULT 1,
--     exp NUMERIC NOT NULL DEFAULT 0
-- );
 CREATE TABLE IF NOT EXISTS tags (
     tag_id SERIAL PRIMARY KEY,
     name VARCHAR(255) NOT NULL,
     location_id BIGINT NOT NULL, -- Typically the Guild ID
     owner_id BIGINT NOT NULL,     -- Discord User ID
     content TEXT,
     attachment_url TEXT,
     tag_type VARCHAR(50) DEFAULT 'default',
     created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
 );
 
 -- 1.a. Create a unique index for LOWER(name), location_id, owner_id
 CREATE UNIQUE INDEX IF NOT EXISTS idx_tags_lower_name_location_owner 
 ON tags (LOWER(name), location_id, owner_id);
 
 -- 2. Create the 'loop_configs' table
 CREATE TABLE IF NOT EXISTS loop_configs (
     guild_id BIGINT PRIMARY KEY,
     channel_id BIGINT,
     enabled BOOLEAN DEFAULT FALSE,
     updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
 );
 
 -- 3. Create the 'borrowed_tags' table
 CREATE TABLE IF NOT EXISTS borrowed_tags (
     borrow_id SERIAL PRIMARY KEY,
     borrower_id BIGINT NOT NULL,       -- Discord User ID who borrowed the tag
     original_tag_id INT NOT NULL,      -- References the original tag
     borrowed_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
     returned_at TIMESTAMP WITH TIME ZONE,
     FOREIGN KEY (original_tag_id) REFERENCES tags(tag_id) ON DELETE CASCADE
 );
 
 -- 3.a. Create a partial unique index to ensure one active borrow per user per tag
 CREATE UNIQUE INDEX IF NOT EXISTS idx_borrowed_tags_unique_active 
 ON borrowed_tags (borrower_id, original_tag_id) 
 WHERE returned_at IS NULL;
 
 -- 4. Create the 'users' table
 CREATE TABLE IF NOT EXISTS users (
     id BIGINT PRIMARY KEY,
     name VARCHAR(255) NOT NULL,
     create_date TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
     level INT NOT NULL DEFAULT 1,
     exp NUMERIC NOT NULL DEFAULT 0
 );
