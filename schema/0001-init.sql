--
-- PostgreSQL database dump
--

-- Dumped from database version 14.20 (Homebrew)
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
    reset boolean DEFAULT false
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
    category text NOT NULL,
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
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.active_flags OWNER TO vyrtuous;

--
-- Name: active_server_voice_mutes; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_server_voice_mutes (
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    expires_in timestamp with time zone,
    member_snowflake bigint NOT NULL,
    reason text,
    updated_at timestamp with time zone DEFAULT now()
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
    reset boolean DEFAULT false
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
    expired boolean DEFAULT false
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
    role_snowflakes bigint[] NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.administrators OWNER TO vyrtuous;

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
    channel_snowflake bigint DEFAULT '-1'::integer,
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    role_snowflake bigint,
    updated_at timestamp with time zone DEFAULT now(),
    CONSTRAINT command_aliases_category_check CHECK ((category = ANY (ARRAY['vegan'::text, 'vmute'::text, 'ban'::text, 'flag'::text, 'tmute'::text, 'role'::text])))
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
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.coordinators OWNER TO vyrtuous;

--
-- Name: developers; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.developers (
    created_at timestamp with time zone DEFAULT now(),
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.developers OWNER TO vyrtuous;

--
-- Name: guild_owners; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.guild_owners (
    created_at timestamp with time zone DEFAULT now(),
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.guild_owners OWNER TO vyrtuous;

--
-- Name: moderation_logs; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.moderation_logs (
    infraction_type text NOT NULL,
    channel_snowflake bigint,
    executor_member_snowflake bigint,
    expires_at timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    is_modification boolean DEFAULT false NOT NULL,
    target_member_snowflake bigint,
    reason text,
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    updated_at timestamp with time zone DEFAULT now() NOT NULL,
    channel_members_voice_count integer DEFAULT 0,
    guild_members_offline_and_online_member_count integer DEFAULT 0,
    guild_members_online_count integer DEFAULT 0,
    guild_members_voice_count integer DEFAULT 0,
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
    updated_at timestamp with time zone DEFAULT now()
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
-- Name: streaming; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.streaming (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    enabled boolean DEFAULT false,
    entry_type text NOT NULL,
    guild_snowflake bigint NOT NULL,
    snowflakes bigint[],
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.streaming OWNER TO vyrtuous;

--
-- Name: sysadmin; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.sysadmin (
    created_at timestamp with time zone DEFAULT now(),
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.sysadmin OWNER TO vyrtuous;

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
-- Name: vegans; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.vegans (
    created_at timestamp with time zone DEFAULT now() NOT NULL,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
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
-- Name: streaming streaming_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.streaming
    ADD CONSTRAINT streaming_pkey PRIMARY KEY (channel_snowflake, entry_type);


--
-- Name: sysadmin sysadmin_member_unique; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.sysadmin
    ADD CONSTRAINT sysadmin_member_unique UNIQUE (member_snowflake);


--
-- Name: temporary_rooms temporary_rooms_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.temporary_rooms
    ADD CONSTRAINT temporary_rooms_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, room_name);


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
-- PostgreSQL database dump complete
--

