DROP TABLE IF EXISTS annotations CASCADE;

DROP TABLE IF EXISTS borrowed_tags CASCADE;

DROP TABLE IF EXISTS citations CASCADE;

DROP TABLE IF EXISTS loop_configs CASCADE;

DROP TABLE IF EXISTS moderation_counts CASCADE;

DROP TABLE IF EXISTS patreon_data CASCADE;

DROP TABLE IF EXISTS pdf_catalog CASCADE;

DROP TABLE IF EXISTS pdfs CASCADE;

DROP TABLE IF EXISTS reference_list CASCADE;

DROP TABLE IF EXISTS reference_tags CASCADE;

DROP TABLE IF EXISTS tags CASCADE;

DROP TABLE IF EXISTS token_usage_logs CASCADE;

DROP TABLE IF EXISTS users CASCADE;

DROP TABLE IF EXISTS faction_members CASCADE;

DROP TABLE IF EXISTS factions CASCADE;

CREATE TABLE public.annotations (
    id integer NOT NULL,
    pdf_id integer,
    user_id bigint NOT NULL,
    page_number integer,
    content text NOT NULL,
    highlighted_text text,
    created_at timestamp without time zone DEFAULT now()
);

CREATE TABLE public.borrowed_tags (
    borrow_id integer NOT NULL,
    borrower_id bigint NOT NULL,
    original_tag_id integer NOT NULL,
    borrowed_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP,
    returned_at timestamp with time zone
);

CREATE TABLE public.citations (
    id integer NOT NULL,
    reference_id integer,
    user_id bigint NOT NULL,
    citation_style text NOT NULL,
    citation_text text NOT NULL,
    created_at timestamp without time zone DEFAULT now()
);

CREATE TABLE public.loop_configs (
    guild_id bigint NOT NULL,
    channel_id bigint,
    enabled boolean DEFAULT false,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE public.pdf_catalog (
    id integer NOT NULL,
    user_id bigint NOT NULL,
    title text NOT NULL,
    file_url text NOT NULL,
    description text,
    tags text[],
    uploaded_at timestamp without time zone DEFAULT now()
);

CREATE TABLE public.tags (
    id integer NOT NULL,
    name varchar(255) NOT NULL,
    location_id bigint NOT NULL,
    owner_id bigint NOT NULL,
    content text,
    attachment_url text,
    tag_type varchar(50) DEFAULT 'default',
    created_at timestamp without time zone DEFAULT now(),
    UNIQUE(name, location_id, owner_id)
);

CREATE TABLE public.users (
    id bigint NOT NULL,
    name varchar(255) NOT NULL,
    create_date timestamp with time zone DEFAULT CURRENT_TIMESTAMP,
    level integer DEFAULT 1 NOT NULL,
    exp numeric DEFAULT 0 NOT NULL
);

CREATE TABLE IF NOT EXISTS public.users (
    id BIGINT PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    create_date TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    level INTEGER DEFAULT 1 NOT NULL,
    exp NUMERIC DEFAULT 0 NOT NULL,
    faction_name VARCHAR(255) DEFAULT NULL
);

-- Create the Factions Table
CREATE TABLE IF NOT EXISTS public.factions (
    name VARCHAR(255) PRIMARY KEY,
    xp NUMERIC DEFAULT 0 NOT NULL,
    level INTEGER DEFAULT 1 NOT NULL
);

-- Create the Faction Members Table (many-to-many relationship)
CREATE TABLE IF NOT EXISTS public.faction_members (
    user_id BIGINT REFERENCES users(id) ON DELETE CASCADE,
    faction_name VARCHAR(255) REFERENCES factions(name) ON DELETE CASCADE,
    PRIMARY KEY (user_id, faction_name)
);
