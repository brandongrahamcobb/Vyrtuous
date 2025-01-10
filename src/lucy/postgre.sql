-- 1. Create the 'tags' table
CREATE TABLE IF NOT EXISTS tags (
    tag_id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    location_id BIGINT NOT NULL, -- Typically the Guild ID
    owner_id BIGINT NOT NULL,     -- Discord User ID
    content TEXT,
    attachment_url TEXT,
    tag_type VARCHAR(50) DEFAULT 'default',
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (LOWER(name), location_id, owner_id)
);

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
    FOREIGN KEY (original_tag_id) REFERENCES tags(tag_id) ON DELETE CASCADE,
    UNIQUE (borrower_id, original_tag_id, returned_at)
);

-- 4. Assuming a 'users' table exists, if not, create it accordingly
CREATE TABLE IF NOT EXISTS users (
    user_id BIGINT PRIMARY KEY,
    username VARCHAR(255) NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);
