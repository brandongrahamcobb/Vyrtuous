-- Drop existing tables
DROP TABLE IF EXISTS annotations CASCADE;
DROP TABLE IF EXISTS pdfs CASCADE;
DROP TABLE IF EXISTS citations CASCADE;
DROP TABLE IF EXISTS reference_list CASCADE;
DROP TABLE IF EXISTS tags CASCADE;
DROP TABLE IF EXISTS loop_configs CASCADE;

-- Create updated schema

CREATE TABLE reference_list (
    id SERIAL PRIMARY KEY,
    user_id BIGINT NOT NULL,
    location_id INT NOT NULL,
    title TEXT NOT NULL,
    authors TEXT[] NOT NULL,
    publication_year INT,
    doi TEXT,
    abstract TEXT,
    tags INT[],
    created_at TIMESTAMP DEFAULT NOW(),
    updated_at TIMESTAMP DEFAULT NOW()
);

CREATE TABLE citations (
    id SERIAL PRIMARY KEY,
    reference_id INT REFERENCES reference_list(id) ON DELETE CASCADE,
    user_id BIGINT NOT NULL,
    citation_style TEXT NOT NULL,
    citation_text TEXT NOT NULL,
    created_at TIMESTAMP DEFAULT NOW()
);

CREATE TABLE pdfs (
    id SERIAL PRIMARY KEY,
    reference_id INT REFERENCES reference_list(id) ON DELETE CASCADE,
    file_url TEXT NOT NULL,
    created_at TIMESTAMP DEFAULT NOW()
);

CREATE TABLE annotations (
    id SERIAL PRIMARY KEY,
    pdf_id INT REFERENCES pdfs(id) ON DELETE CASCADE,
    user_id BIGINT NOT NULL,
    page_number INT,
    content TEXT NOT NULL,
    highlighted_text TEXT,
    created_at TIMESTAMP DEFAULT NOW()
);

-- Existing Tag Table
CREATE TABLE IF NOT EXISTS tags (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    location_id BIGINT NOT NULL,
    content TEXT,
    attachment_url TEXT,
    owner_id BIGINT NOT NULL,
    tag_type VARCHAR(50) DEFAULT 'default',
    UNIQUE(name, location_id, owner_id)
);

-- Loop Config Table (Assuming it exists from tag.py)
CREATE TABLE IF NOT EXISTS loop_configs (
    guild_id BIGINT PRIMARY KEY,
    channel_id BIGINT,
    enabled BOOLEAN DEFAULT TRUE
);
