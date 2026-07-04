--
-- PostgreSQL database dump
--

\restrict FNtfa75ifb0GgUBMS2Z7O80ZHWVbi7NjTGhKLKQcnDi8YPtKzMFyG1uxWQ3oY5q

-- Dumped from database version 18.4 (Debian 18.4-1.pgdg13+1)
-- Dumped by pg_dump version 18.4

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
-- Name: active_automute_channels; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_automute_channels (
    channel_snowflake bigint CONSTRAINT active_stages_channel_snowflake_not_null NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint CONSTRAINT active_stages_guild_snowflake_not_null NOT NULL,
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
    reason text,
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
    category text CONSTRAINT active_caps_moderation_type_not_null NOT NULL,
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
    reason text,
    updated_at timestamp with time zone DEFAULT now()
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
    updated_at timestamp with time zone DEFAULT now()
);


ALTER TABLE public.active_server_voice_mutes OWNER TO vyrtuous;

--
-- Name: active_text_mutes; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_text_mutes (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text,
    updated_at timestamp with time zone DEFAULT now(),
    last_muted timestamp with time zone DEFAULT now(),
    reset boolean DEFAULT false
);


ALTER TABLE public.active_text_mutes OWNER TO vyrtuous;

--
-- Name: active_voice_mutes; Type: TABLE; Schema: public; Owner: vyrtuous
--

CREATE TABLE public.active_voice_mutes (
    channel_snowflake bigint NOT NULL,
    created_at timestamp with time zone DEFAULT now(),
    expires_in timestamp with time zone,
    guild_snowflake bigint NOT NULL,
    member_snowflake bigint NOT NULL,
    reason text,
    target text,
    updated_at timestamp with time zone DEFAULT now()
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
    target_highest_role text,
    target text DEFAULT 'user'::text,
    role_snowflake bigint
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
-- Name: streaming id; Type: DEFAULT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.streaming ALTER COLUMN id SET DEFAULT nextval('public.history_id_seq'::regclass);


--
-- Data for Name: active_automute_channels; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.active_automute_channels (channel_snowflake, created_at, expires_in, guild_snowflake, updated_at) FROM stdin;
\.


--
-- Data for Name: active_bans; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.active_bans (channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, updated_at, last_kicked, reset, blacklisted) FROM stdin;
\.


--
-- Data for Name: active_caps; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.active_caps (channel_snowflake, created_at, duration_seconds, guild_snowflake, category, updated_at) FROM stdin;
10000000000000011	2026-07-02 13:44:09.414674+00	28800	10000000000000500	ban	2026-07-02 13:44:09.414675+00
\.


--
-- Data for Name: active_flags; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.active_flags (channel_snowflake, created_at, guild_snowflake, member_snowflake, reason, updated_at) FROM stdin;
\.


--
-- Data for Name: active_members; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.active_members (created_at, display_name, guild_snowflake, last_active, member_snowflake, updated_at) FROM stdin;
\.


--
-- Data for Name: active_server_voice_mutes; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.active_server_voice_mutes (created_at, guild_snowflake, expires_in, member_snowflake, reason, updated_at) FROM stdin;
\.


--
-- Data for Name: active_text_mutes; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.active_text_mutes (channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, updated_at, last_muted, reset) FROM stdin;
\.


--
-- Data for Name: active_voice_mutes; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.active_voice_mutes (channel_snowflake, created_at, expires_in, guild_snowflake, member_snowflake, reason, target, updated_at) FROM stdin;
\.


--
-- Data for Name: administrator_roles; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.administrator_roles (created_at, guild_snowflake, role_snowflake, updated_at) FROM stdin;
\.


--
-- Data for Name: administrators; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.administrators (created_at, guild_snowflake, member_snowflake, role_snowflakes, updated_at) FROM stdin;
\.


--
-- Data for Name: bug_tracking; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.bug_tracking (channel_snowflake, created_at, member_snowflakes, guild_snowflake, id, message_snowflake, notes, resolved, updated_at, expired) FROM stdin;
\.


--
-- Data for Name: command_aliases; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.command_aliases (category, alias_name, channel_snowflake, created_at, guild_snowflake, role_snowflake, updated_at) FROM stdin;
ban	vban	1347284829305311273	2026-07-01 15:32:33.502862+00	1347284827350630591	\N	2026-07-01 15:32:33.502864+00
vmute	testmute	10000000000000011	2026-07-04 00:54:54.796775+00	10000000000000500	\N	2026-07-04 00:54:54.796777+00
flag	testflag	10000000000000011	2026-07-04 00:54:55.626229+00	10000000000000500	\N	2026-07-04 00:54:55.626232+00
vegan	testvegan	10000000000000011	2026-07-04 00:54:56.419884+00	10000000000000500	\N	2026-07-04 00:54:56.419886+00
tmute	testtmute	10000000000000011	2026-07-04 00:54:57.235411+00	10000000000000500	10000000000000200	2026-07-04 00:54:57.235413+00
ban	testban	10000000000000011	2026-07-04 00:54:58.075006+00	10000000000000500	10000000000000200	2026-07-04 00:54:58.075009+00
role	testrole	10000000000000011	2026-07-04 00:54:58.956548+00	10000000000000500	10000000000000200	2026-07-04 00:54:58.95655+00
\.


--
-- Data for Name: coordinators; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.coordinators (channel_snowflake, created_at, guild_snowflake, member_snowflake, updated_at) FROM stdin;
\.


--
-- Data for Name: developers; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.developers (created_at, member_snowflake, updated_at) FROM stdin;
\.


--
-- Data for Name: moderation_logs; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.moderation_logs (identifier, channel_snowflake, author_snowflake, expires_at, guild_snowflake, is_modification, target_snowflake, reason, created_at, updated_at, current_channel_members, total_guild_members, online_members, total_voice_members, executor_highest_role, target_highest_role, target, role_snowflake) FROM stdin;
ban	\N	154749533429956608	2026-07-02 00:03:07.253796+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-01 16:03:07.253894+00	2026-07-01 16:03:07.253894+00	0	0	0	0	Sysadmin	Everyone	user	\N
ban	\N	154749533429956608	2026-07-02 21:44:12.165139+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 13:44:12.165268+00	2026-07-02 13:44:12.165268+00	0	0	0	0	Sysadmin	Everyone	user	\N
unban	\N	154749533429956608	2026-07-02 13:46:26.274123+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 13:46:26.27417+00	2026-07-02 13:46:26.27417+00	0	0	0	0	Sysadmin	Everyone	user	\N
ban	\N	154749533429956608	2026-07-02 21:48:21.96745+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 13:48:21.96757+00	2026-07-02 13:48:21.96757+00	0	0	0	0	Sysadmin	Everyone	user	\N
unban	\N	154749533429956608	2026-07-02 13:49:26.825286+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 13:49:26.825329+00	2026-07-02 13:49:26.825329+00	0	0	0	0	Sysadmin	Everyone	user	\N
ban	\N	154749533429956608	2026-07-02 21:50:57.710692+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 13:50:57.710831+00	2026-07-02 13:50:57.710831+00	0	0	0	0	Sysadmin	Everyone	user	\N
unban	\N	154749533429956608	2026-07-02 13:51:00.651436+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 13:51:00.651536+00	2026-07-02 13:51:00.651536+00	0	0	0	0	Sysadmin	Everyone	user	\N
ban	\N	154749533429956608	2026-07-02 15:16:42.704248+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 15:08:42.704336+00	2026-07-02 15:08:42.704336+00	0	0	0	0	Sysadmin	Everyone	user	\N
unban	\N	154749533429956608	2026-07-02 15:08:51.037376+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 15:08:51.037424+00	2026-07-02 15:08:51.037424+00	0	0	0	0	Sysadmin	Everyone	user	\N
ban	\N	154749533429956608	2026-07-02 15:16:55.669726+00	1347284827350630591	f	1454139068852994090	test	2026-07-02 15:08:55.669785+00	2026-07-02 15:08:55.669785+00	0	0	0	0	Sysadmin	Everyone	user	\N
unban	\N	154749533429956608	2026-07-02 15:09:00.783414+00	1347284827350630591	f	1454139068852994090	No reason provided.	2026-07-02 15:09:00.783478+00	2026-07-02 15:09:00.783478+00	0	0	0	0	Sysadmin	Everyone	user	\N
unvmute	\N	\N	2026-07-03 18:55:05.638109+00	1347284827350630591	f	154749533429956608	No reason provided.	2026-07-03 18:55:05.638497+00	2026-07-03 18:55:05.638497+00	1	0	0	0	Role undetermined	Role undetermined	user	\N
vmute	\N	\N	2026-07-03 19:55:08.458649+00	1347284827350630591	f	154749533429956608	Right-click muted	2026-07-03 18:55:08.459077+00	2026-07-03 18:55:08.459077+00	1	0	0	0	Role undetermined	Role undetermined	user	\N
unvmute	\N	\N	2026-07-03 18:55:10.658727+00	1347284827350630591	f	154749533429956608	No reason provided.	2026-07-03 18:55:10.658917+00	2026-07-03 18:55:10.658917+00	1	0	0	0	Role undetermined	Role undetermined	user	\N
vmute	\N	\N	2026-07-03 19:55:13.845321+00	1347284827350630591	f	154749533429956608	Right-click muted	2026-07-03 18:55:13.845541+00	2026-07-03 18:55:13.845541+00	1	0	0	0	Role undetermined	Role undetermined	user	\N
vmute	\N	\N	2026-07-03 19:55:14.054887+00	1347284827350630591	f	154749533429956608	Right-click muted	2026-07-03 18:55:14.055031+00	2026-07-03 18:55:14.055031+00	1	0	0	0	Role undetermined	Role undetermined	user	\N
\.


--
-- Data for Name: moderators; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.moderators (channel_snowflake, created_at, guild_snowflake, member_snowflake, updated_at) FROM stdin;
\.


--
-- Data for Name: streaming; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.streaming (target_channel_snowflake, created_at, guild_snowflake, id, updated_at, source_channel_snowflake) FROM stdin;
\.


--
-- Data for Name: temporary_rooms; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.temporary_rooms (channel_snowflake, created_at, guild_snowflake, room_name, updated_at) FROM stdin;
\.


--
-- Data for Name: uploads; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.uploads (command_name, file_bytes, filename, tag, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: vegans; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.vegans (created_at, guild_snowflake, member_snowflake, updated_at) FROM stdin;
\.


--
-- Data for Name: video_rooms; Type: TABLE DATA; Schema: public; Owner: vyrtuous
--

COPY public.video_rooms (channel_snowflake, created_at, guild_snowflake, updated_at) FROM stdin;
\.


--
-- Name: history_id_seq; Type: SEQUENCE SET; Schema: public; Owner: vyrtuous
--

SELECT pg_catalog.setval('public.history_id_seq', 1, false);


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
-- Name: active_automute_channels active_stages_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.active_automute_channels
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
-- Name: moderators moderators_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.moderators
    ADD CONSTRAINT moderators_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, member_snowflake);


--
-- Name: temporary_rooms temporary_rooms_pkey; Type: CONSTRAINT; Schema: public; Owner: vyrtuous
--

ALTER TABLE ONLY public.temporary_rooms
    ADD CONSTRAINT temporary_rooms_pkey PRIMARY KEY (channel_snowflake, guild_snowflake, room_name);


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

\unrestrict FNtfa75ifb0GgUBMS2Z7O80ZHWVbi7NjTGhKLKQcnDi8YPtKzMFyG1uxWQ3oY5q

