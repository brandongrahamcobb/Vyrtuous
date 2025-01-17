```sql
-- template_database.sql
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

-- References Table
CREATE TABLE IF NOT EXISTS references (
    id SERIAL PRIMARY KEY,
    user_id BIGINT NOT NULL,
    location_id BIGINT NOT NULL,
    title VARCHAR(1024) NOT NULL,
    authors TEXT[],
    publication_year INT,
    doi VARCHAR(255),
    abstract TEXT,
    tags INT[] REFERENCES tags(id) ON DELETE SET NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Index for faster search on title and authors
CREATE INDEX IF NOT EXISTS idx_references_title ON references USING GIN (to_tsvector('english', title));
CREATE INDEX IF NOT EXISTS idx_references_authors ON references USING GIN (authors);

-- PDFs Table
CREATE TABLE IF NOT EXISTS pdfs (
    id SERIAL PRIMARY KEY,
    reference_id INT REFERENCES references(id) ON DELETE CASCADE,
    file_url TEXT NOT NULL,
    uploaded_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Annotations Table
CREATE TABLE IF NOT EXISTS annotations (
    id SERIAL PRIMARY KEY,
    pdf_id INT REFERENCES pdfs(id) ON DELETE CASCADE,
    user_id BIGINT NOT NULL,
    page_number INT,
    content TEXT,
    highlighted_text TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Citations Table
CREATE TABLE IF NOT EXISTS citations (
    id SERIAL PRIMARY KEY,
    reference_id INT REFERENCES references(id) ON DELETE CASCADE,
    user_id BIGINT NOT NULL,
    citation_style VARCHAR(50) NOT NULL,
    citation_text TEXT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```