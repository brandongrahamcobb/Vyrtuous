--
-- PostgreSQL database dump
--

-- Dumped from database version 14.23 (Homebrew)
-- Dumped by pg_dump version 14.18 (Homebrew)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: postgres_fdw; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS postgres_fdw WITH SCHEMA public;


--
-- Name: EXTENSION postgres_fdw; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION postgres_fdw IS 'foreign-data wrapper for remote PostgreSQL servers';


--
-- Name: set_expired(); Type: FUNCTION; Schema: public; Owner: vyrtuous
--

CREATE FUNCTION public.set_expired() RETURNS trigger
    LANGUAGE plpgsql
    AS $$
BEGIN
    NEW.expired := (NEW.expires_in IS NOT NULL AND NEW.expires_in < NOW());
    RETURN NEW;
END;
$$;


ALTER FUNCTION public.set_expired() OWNER TO vyrtuous;

--
-- Name: tmpdb; Type: SERVER; Schema: -; Owner: vyrtuous
--

CREATE SERVER tmpdb FOREIGN DATA WRAPPER postgres_fdw OPTIONS (
    dbname 'vyrtuous_tmp'
);


ALTER SERVER tmpdb OWNER TO vyrtuous;

--
-- Name: USER MAPPING vyrtuous SERVER tmpdb; Type: USER MAPPING; Schema: -; Owner: vyrtuous
--

CREATE USER MAPPING FOR vyrtuous SERVER tmpdb OPTIONS (
    "user" 'vyrtuous'
);


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: active_automute_channels; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_automute_channels (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.active_automute_channels OWNER TO vyrtuous;

--
-- Name: active_bans; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_bans (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    last_kicked timestamp with time zone DEFAULT now(),
    reset boolean DEFAULT false,
    blacklisted boolean
);


ALTER TABLE public.active_bans OWNER TO vyrtuous;

--
-- Name: active_caps; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_caps (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    duration_seconds integer NOT NULL,
    guild_snowflake bigint NOT NULL,
    category text NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    CONSTRAINT active_caps_moderation_type_check CHECK ((category = ANY (ARRAY['ban'::text, 'vmute'::text, 'tmute'::text])))
);


ALTER TABLE public.active_caps OWNER TO vyrtuous;

--
-- Name: active_flags; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_flags (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.active_flags OWNER TO vyrtuous;

--
-- Name: active_members; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_members (
    created_at timestamp with time zone DEFAULT now(),
    display_name text,
    last_active timestamp with time zone DEFAULT now(),
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.active_members OWNER TO vyrtuous;

--
-- Name: active_text_mutes; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_text_mutes (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    last_muted timestamp with time zone DEFAULT now(),
    reset boolean DEFAULT false
);


ALTER TABLE public.active_text_mutes OWNER TO vyrtuous;

--
-- Name: active_video_only_channels; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_video_only_channels (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone
);


ALTER TABLE public.active_video_only_channels OWNER TO vyrtuous;

--
-- Name: active_voice_mutes; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_voice_mutes (
    channel_snowflake bigint,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text NOT NULL,
    target text NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    id integer NOT NULL
);


ALTER TABLE public.active_voice_mutes OWNER TO vyrtuous;

--
-- Name: active_voice_mutes_id_seq; Type: SEQUENCE; Schema: public; Owner: vyrtuous
--

CREATE SEQUENCE public.active_voice_mutes_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.active_voice_mutes_id_seq OWNER TO vyrtuous;

--
-- Name: active_voice_mutes_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: vyrtuous
--

ALTER SEQUENCE public.active_voice_mutes_id_seq OWNED BY public.active_voice_mutes.id;


--
-- Name: autoassign_roles; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.autoassign_roles (
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    group_alias text NOT NULL,
    channel_snowflake bigint
);


ALTER TABLE public.autoassign_roles OWNER TO vyrtuous;

--
-- Name: bug_tracking; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.bug_tracking (
    channel_snowflake bigint,
    created_at timestamp with time zone DEFAULT now(),
    member_snowflakes bigint[],
    guild_snowflake bigint NOT NULL,
    id uuid NOT NULL,
    message_snowflake bigint NOT NULL,
    notes text,
    resolved boolean DEFAULT false NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    expired boolean DEFAULT false
);


ALTER TABLE public.bug_tracking OWNER TO vyrtuous;

--
-- Name: command_aliases; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.command_aliases (
    category text NOT NULL,
    alias_name text NOT NULL,
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint,
    updated_at timestamp with time zone DEFAULT now(),
    CONSTRAINT command_aliases_category_check CHECK ((category = ANY (ARRAY['vmute'::text, 'ban'::text, 'flag'::text, 'tmute'::text, 'role'::text])))
);


ALTER TABLE public.command_aliases OWNER TO vyrtuous;

--
-- Name: moderation_logs; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.moderation_logs (
    identifier text NOT NULL,
    channel_snowflake bigint,
    author_snowflake bigint,
    expires_at timestamp with time zone,
    guild_snowflake bigint,
    is_modification boolean DEFAULT false NOT NULL,
    target_snowflake bigint,
    reason text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL,
    current_channel_members integer DEFAULT 0,
    total_guild_members integer DEFAULT 0,
    online_members integer DEFAULT 0,
    total_voice_members integer DEFAULT 0,
    executor_highest_role text,
    target_highest_role text,
    target text DEFAULT 'user'::text,
    role_snowflake bigint,
    id integer NOT NULL
);


ALTER TABLE public.moderation_logs OWNER TO vyrtuous;

--
-- Name: moderation_logs_id_seq; Type: SEQUENCE; Schema: public; Owner: vyrtuous
--

CREATE SEQUENCE public.moderation_logs_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.moderation_logs_id_seq OWNER TO vyrtuous;

--
-- Name: moderation_logs_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: vyrtuous
--

ALTER SEQUENCE public.moderation_logs_id_seq OWNED BY public.moderation_logs.id;


--
-- Name: permission_entries; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.permission_entries (
    channel_snowflake bigint,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint,
    group_alias text NOT NULL,
    member_snowflake bigint NOT NULL,
    role_snowflakes bigint[] NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.permission_entries OWNER TO vyrtuous;

--
-- Name: roles; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.roles (
    created_at timestamp with time zone DEFAULT now(),
    channel_snowflake bigint NOT NULL,
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.roles OWNER TO vyrtuous;

--
-- Name: streaming; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.streaming (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    source_guild_snowflake bigint,
    id bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    source_channel_snowflake bigint,
    guild_snowflake bigint NOT NULL
);


ALTER TABLE public.streaming OWNER TO vyrtuous;

--
-- Name: streaming_id_seq; Type: SEQUENCE; Schema: public; Owner: vyrtuous
--

CREATE SEQUENCE public.streaming_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.streaming_id_seq OWNER TO vyrtuous;

--
-- Name: streaming_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: vyrtuous
--

ALTER SEQUENCE public.streaming_id_seq OWNED BY public.streaming.id;


--
-- Name: uploads; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.uploads (
    command_name text NOT NULL,
    file_bytes bytea NOT NULL,
    filename text NOT NULL,
    tag text NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.uploads OWNER TO vyrtuous;

--
-- Name: vegans; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.vegans (
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    notes text NOT NULL
);


ALTER TABLE public.vegans OWNER TO vyrtuous;

--
-- Name: active_voice_mutes id; Type: DEFAULT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_voice_mutes ALTER COLUMN id SET DEFAULT nextval('public.active_voice_mutes_id_seq'::regclass);


--
-- Name: moderation_logs id; Type: DEFAULT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.moderation_logs ALTER COLUMN id SET DEFAULT nextval('public.moderation_logs_id_seq'::regclass);


--
-- Name: streaming id; Type: DEFAULT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.streaming ALTER COLUMN id SET DEFAULT nextval('public.streaming_id_seq'::regclass);


--
-- Name: active_automute_channels active_automute_channels_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_automute_channels
    ADD CONSTRAINT active_automute_channels_pkey PRIMARY KEY (guild_snowflake, channel_snowflake);


--
-- Name: active_bans active_bans_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_bans
    ADD CONSTRAINT active_bans_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake);


--
-- Name: active_caps active_caps_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_caps
    ADD CONSTRAINT active_caps_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, category);


--
-- Name: active_flags active_flags_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_flags
    ADD CONSTRAINT active_flags_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake);


--
-- Name: active_members active_members_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_members
    ADD CONSTRAINT active_members_pkey PRIMARY KEY (member_snowflake);


--
-- Name: active_text_mutes active_text_mutes_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_text_mutes
    ADD CONSTRAINT active_text_mutes_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake);


--
-- Name: active_video_only_channels active_video_only_channels_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_video_only_channels
    ADD CONSTRAINT active_video_only_channels_pkey PRIMARY KEY (channel_snowflake, guild_snowflake);

--
-- Name: active_video_only_channels active_video_only_channels_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_voice_mutes
    ADD CONSTRAINT active_voice_mutes_pkey PRIMARY KEY (id);

--
-- Name: autoassign_roles autoassign_roles_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.autoassign_roles
    ADD CONSTRAINT autoassign_roles_pkey PRIMARY KEY (guild_snowflake, role_snowflake);


--
-- Name: bug_tracking bug_tracking_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.bug_tracking
    ADD CONSTRAINT bug_tracking_pkey PRIMARY KEY (id);


--
-- Name: command_aliases command_aliases_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.command_aliases
    ADD CONSTRAINT command_aliases_pkey PRIMARY KEY (alias_name, category, guild_snowflake);


--
-- Name: moderation_logs moderation_logs_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.moderation_logs
    ADD CONSTRAINT moderation_logs_pkey PRIMARY KEY (id);


--
-- Name: permission_entries permission_entries_unique; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.permission_entries
    ADD CONSTRAINT permission_entries_unique UNIQUE (member_snowflake, group_alias, guild_snowflake, channel_snowflake);


--
-- Name: streaming unique_target_source; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.streaming
    ADD CONSTRAINT unique_target_source UNIQUE (channel_snowflake, source_channel_snowflake);


--
-- Name: uploads uploads_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.uploads
    ADD CONSTRAINT uploads_pkey PRIMARY KEY (command_name, tag);


--
-- Name: vegans vegans_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.vegans
    ADD CONSTRAINT vegans_pkey PRIMARY KEY (guild_snowflake, member_snowflake);


--
-- PostgreSQL database dump complete
--

