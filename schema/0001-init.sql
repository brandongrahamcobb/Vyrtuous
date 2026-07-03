--
-- PostgreSQL database dump
--

\restrict taP7PpPY8bZmusHt1Qnm4cvtklEmjha0mSSa8TPstkM4xE0WDzTbTZnqv69qh5t

-- Dumped from database version 18.0 (Debian 18.0-1.pgdg13+3)
-- Dumped by pg_dump version 18.0 (Debian 18.0-1.pgdg13+3)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET transaction_timeout = 0;
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
-- Name: active_bans; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_bans (
    channel_snowflake bigint DEFAULT '-1'::integer NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text,
    updated_at timestamp with time zone DEFAULT now(),
    expired boolean DEFAULT false,
    last_kicked timestamp with time zone DEFAULT now(),
    reset boolean DEFAULT false,
    display_name text,
    blacklisted boolean
);


ALTER TABLE public.active_bans OWNER TO vyrtuous;

--
-- Name: active_caps; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_caps (
    channel_snowflake bigint DEFAULT '-1'::integer NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    duration_seconds integer NOT NULL,
    guild_snowflake bigint NOT NULL,
    category text CONSTRAINT active_caps_moderation_type_not_null NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    CONSTRAINT active_caps_moderation_type_check CHECK ((category = ANY (ARRAY['ban'::text, 'vmute'::text, 'tmute'::text])))
);


ALTER TABLE public.active_caps OWNER TO vyrtuous;

--
-- Name: active_flags; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_flags (
    channel_snowflake bigint DEFAULT '-1'::integer NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text
);


ALTER TABLE public.active_flags OWNER TO vyrtuous;

--
-- Name: active_members; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_members (
    created_at timestamp with time zone DEFAULT now(),
    display_name text,
    guild_snowflake bigint NOT NULL,
    last_active timestamp with time zone DEFAULT now(),
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.active_members OWNER TO vyrtuous;

--
-- Name: active_server_voice_mutes; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_server_voice_mutes (
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    expires_in timestamp with time zone,
    member_snowflake bigint NOT NULL,
    reason text,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text
);


ALTER TABLE public.active_server_voice_mutes OWNER TO vyrtuous;

--
-- Name: active_stages; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_stages (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    expired boolean DEFAULT false
);


ALTER TABLE public.active_stages OWNER TO vyrtuous;

--
-- Name: active_text_mutes; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_text_mutes (
    channel_snowflake bigint DEFAULT '-1'::integer NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text,
    updated_at timestamp with time zone DEFAULT now(),
    expired boolean DEFAULT false,
    role_snowflake bigint,
    last_muted timestamp with time zone DEFAULT now(),
    reset boolean DEFAULT false,
    display_name text
);


ALTER TABLE public.active_text_mutes OWNER TO vyrtuous;

--
-- Name: active_voice_mutes; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_voice_mutes (
    channel_snowflake bigint DEFAULT '-1'::integer NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text,
    target text,
    updated_at timestamp with time zone DEFAULT now(),
    expired boolean DEFAULT false,
    display_name text
);


ALTER TABLE public.active_voice_mutes OWNER TO vyrtuous;

--
-- Name: administrator_roles; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.administrator_roles (
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.administrator_roles OWNER TO vyrtuous;

--
-- Name: administrators; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.administrators (
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    role_snowflakes bigint[] CONSTRAINT administrators_role_snowflake_not_null NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text
);


ALTER TABLE public.administrators OWNER TO vyrtuous;

--
-- Name: ban_roles; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.ban_roles (
    created_at timestamp with time zone DEFAULT now(),
    channel_snowflake bigint NOT NULL,
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.ban_roles OWNER TO vyrtuous;

--
-- Name: bug_tracking; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.bug_tracking (
    channel_snowflake bigint,
    created_at timestamp with time zone DEFAULT now(),
    member_snowflakes bigint[],
    guild_snowflake bigint CONSTRAINT developer_logs_guild_snowflake_not_null NOT NULL,
    id uuid CONSTRAINT developer_logs_id_not_null NOT NULL,
    message_snowflake bigint CONSTRAINT developer_logs_message_snowflake_not_null NOT NULL,
    notes text,
    resolved boolean DEFAULT false CONSTRAINT developer_logs_resolved_not_null NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    expired boolean DEFAULT false
);


ALTER TABLE public.bug_tracking OWNER TO vyrtuous;

--
-- Name: command_aliases; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.command_aliases (
    category text CONSTRAINT command_aliases_alias_type_not_null NOT NULL,
    alias_name text NOT NULL,
    channel_snowflake bigint DEFAULT '-1'::integer,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint,
    updated_at timestamp with time zone DEFAULT now(),
    CONSTRAINT command_aliases_category_check CHECK ((category = ANY (ARRAY['vegan'::text, 'vmute'::text, 'ban'::text, 'flag'::text, 'tmute'::text, 'role'::text, 'hide'::text])))
);


ALTER TABLE public.command_aliases OWNER TO vyrtuous;

--
-- Name: coordinators; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.coordinators (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text
);


ALTER TABLE public.coordinators OWNER TO vyrtuous;

--
-- Name: developers; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.developers (
    created_at timestamp with time zone DEFAULT now(),
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text
);


ALTER TABLE public.developers OWNER TO vyrtuous;

--
-- Name: guild_owners; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.guild_owners (
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text
);


ALTER TABLE public.guild_owners OWNER TO vyrtuous;

--
-- Name: hide_roles; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.hide_roles (
    created_at timestamp with time zone DEFAULT now(),
    channel_snowflake bigint NOT NULL,
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.hide_roles OWNER TO vyrtuous;

--
-- Name: streaming; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.streaming (
    target_channel_snowflake bigint CONSTRAINT history_channel_snowflake_not_null NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint CONSTRAINT history_guild_snowflake_not_null NOT NULL,
    id bigint CONSTRAINT history_id_not_null NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    source_channel_snowflake bigint
);


ALTER TABLE public.streaming OWNER TO vyrtuous;

--
-- Name: history_id_seq; Type: SEQUENCE; Schema: public; Owner: vyrtuous
--

CREATE SEQUENCE public.history_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.history_id_seq OWNER TO vyrtuous;

--
-- Name: history_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: vyrtuous
--

ALTER SEQUENCE public.history_id_seq OWNED BY public.streaming.id;


--
-- Name: moderation_logs; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.moderation_logs (
    identifier text CONSTRAINT moderation_logs_action_type_not_null NOT NULL,
    channel_snowflake bigint,
    author_snowflake bigint,
    expires_at timestamp with time zone,
    guild_snowflake bigint NOT NULL,
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
    target_highest_role text
);


ALTER TABLE public.moderation_logs OWNER TO vyrtuous;

--
-- Name: moderators; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.moderators (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text
);


ALTER TABLE public.moderators OWNER TO vyrtuous;

--
-- Name: roles; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.roles (
    created_at timestamp with time zone DEFAULT now(),
    channel_snowflake bigint NOT NULL,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    role_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.roles OWNER TO vyrtuous;

--
-- Name: sysadmin; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.sysadmin (
    id boolean DEFAULT true NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text,
    CONSTRAINT sysadmin_id_check CHECK (id)
);


ALTER TABLE public.sysadmin OWNER TO vyrtuous;

--
-- Name: temporary_blacklist; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.temporary_blacklist (
    member_snowflake bigint NOT NULL,
    member_display_name text
);


ALTER TABLE public.temporary_blacklist OWNER TO vyrtuous;

--
-- Name: temporary_rooms; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.temporary_rooms (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    room_name text NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.temporary_rooms OWNER TO vyrtuous;

--
-- Name: text_mute_roles; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.text_mute_roles (
    created_at timestamp with time zone DEFAULT now(),
    channel_snowflake bigint NOT NULL,
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.text_mute_roles OWNER TO vyrtuous;

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
-- Name: users; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.users (
    discord_snowflake bigint CONSTRAINT users_discord_snowflake_not_null1 NOT NULL,
    moderator_channel_ids bigint[],
    coordinator_channel_ids bigint[],
    developer_guild_ids bigint[],
    server_mute_guild_ids bigint[],
    administrator_guild_ids bigint[],
    updated_at timestamp with time zone DEFAULT now(),
    created_at timestamp with time zone DEFAULT now(),
    administrator_role_ids bigint[]
);


ALTER TABLE public.users OWNER TO vyrtuous;

--
-- Name: vegans; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.vegans (
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now(),
    display_name text
);


ALTER TABLE public.vegans OWNER TO vyrtuous;

--
-- Name: video_rooms; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.video_rooms (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.video_rooms OWNER TO vyrtuous;

--
-- Name: streaming id; Type: DEFAULT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.streaming ALTER COLUMN id SET DEFAULT nextval('public.history_id_seq'::regclass);


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
-- Name: active_server_voice_mutes active_server_voice_mutes_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_server_voice_mutes
    ADD CONSTRAINT active_server_voice_mutes_pkey PRIMARY KEY (guild_snowflake, member_snowflake);


--
-- Name: active_stages active_stages_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_stages
    ADD CONSTRAINT active_stages_pkey PRIMARY KEY (guild_snowflake, channel_snowflake);


--
-- Name: active_text_mutes active_text_mutes_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_text_mutes
    ADD CONSTRAINT active_text_mutes_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake);


--
-- Name: active_voice_mutes active_voice_mutes_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_voice_mutes
    ADD CONSTRAINT active_voice_mutes_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake);


--
-- Name: administrator_roles administrator_roles_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.administrator_roles
    ADD CONSTRAINT administrator_roles_pkey PRIMARY KEY (guild_snowflake, role_snowflake);


--
-- Name: administrators administrators_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.administrators
    ADD CONSTRAINT administrators_pkey PRIMARY KEY (guild_snowflake, member_snowflake);


--
-- Name: ban_roles ban_roles_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.ban_roles
    ADD CONSTRAINT ban_roles_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, role_snowflake);


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
-- Name: coordinators coordinators_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.coordinators
    ADD CONSTRAINT coordinators_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake);


--
-- Name: developers developers_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.developers
    ADD CONSTRAINT developers_pkey PRIMARY KEY (member_snowflake);


--
-- Name: guild_owners guild_owners_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.guild_owners
    ADD CONSTRAINT guild_owners_pkey PRIMARY KEY (guild_snowflake, member_snowflake);


--
-- Name: hide_roles hide_roles_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.hide_roles
    ADD CONSTRAINT hide_roles_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, role_snowflake);


--
-- Name: moderators moderators_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.moderators
    ADD CONSTRAINT moderators_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake);


--
-- Name: roles roles_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.roles
    ADD CONSTRAINT roles_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake, role_snowflake);


--
-- Name: sysadmin sysadmin_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.sysadmin
    ADD CONSTRAINT sysadmin_pkey PRIMARY KEY (id);


--
-- Name: temporary_blacklist temporary_blacklist_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.temporary_blacklist
    ADD CONSTRAINT temporary_blacklist_pkey PRIMARY KEY (member_snowflake);


--
-- Name: temporary_rooms temporary_rooms_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.temporary_rooms
    ADD CONSTRAINT temporary_rooms_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, room_name);


--
-- Name: text_mute_roles text_mute_roles_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.text_mute_roles
    ADD CONSTRAINT text_mute_roles_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, role_snowflake);


--
-- Name: streaming unique_target_source; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.streaming
    ADD CONSTRAINT unique_target_source UNIQUE (target_channel_snowflake, source_channel_snowflake);


--
-- Name: uploads uploads_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.uploads
    ADD CONSTRAINT uploads_pkey PRIMARY KEY (command_name, tag);


--
-- Name: users users_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_pkey PRIMARY KEY (discord_snowflake);


--
-- Name: vegans vegans_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.vegans
    ADD CONSTRAINT vegans_pkey PRIMARY KEY (guild_snowflake, member_snowflake);


--
-- Name: video_rooms video_rooms_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.video_rooms
    ADD CONSTRAINT video_rooms_pkey PRIMARY KEY (channel_snowflake, guild_snowflake);


--
-- Name: active_bans set_expired_active_bans; Type: TRIGGER; Schema: public; Owner: vyrtuous
--

CREATE TRIGGER set_expired_active_bans BEFORE INSERT OR UPDATE ON public.active_bans FOR EACH ROW EXECUTE FUNCTION public.set_expired();


--
-- Name: active_stages set_expired_active_stages; Type: TRIGGER; Schema: public; Owner: vyrtuous
--

CREATE TRIGGER set_expired_active_stages BEFORE INSERT OR UPDATE ON public.active_stages FOR EACH ROW EXECUTE FUNCTION public.set_expired();


--
-- Name: active_text_mutes set_expired_active_text_mutes; Type: TRIGGER; Schema: public; Owner: vyrtuous
--

CREATE TRIGGER set_expired_active_text_mutes BEFORE INSERT OR UPDATE ON public.active_text_mutes FOR EACH ROW EXECUTE FUNCTION public.set_expired();


--
-- Name: active_voice_mutes set_expired_active_voice_mutes; Type: TRIGGER; Schema: public; Owner: vyrtuous
--

CREATE TRIGGER set_expired_active_voice_mutes BEFORE INSERT OR UPDATE ON public.active_voice_mutes FOR EACH ROW EXECUTE FUNCTION public.set_expired();


--
-- PostgreSQL database dump complete
--

\unrestrict taP7PpPY8bZmusHt1Qnm4cvtklEmjha0mSSa8TPstkM4xE0WDzTbTZnqv69qh5t

