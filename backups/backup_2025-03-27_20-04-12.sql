--
-- PostgreSQL database dump
--

-- Dumped from database version 17.2
-- Dumped by pg_dump version 17.2

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

SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: annotations; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.annotations (
    id integer NOT NULL,
    pdf_id integer,
    user_id bigint NOT NULL,
    page_number integer,
    content text NOT NULL,
    highlighted_text text,
    created_at timestamp without time zone DEFAULT now()
);


ALTER TABLE public.annotations OWNER TO postgres;

--
-- Name: annotations_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.annotations_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.annotations_id_seq OWNER TO postgres;

--
-- Name: annotations_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.annotations_id_seq OWNED BY public.annotations.id;


--
-- Name: borrowed_tags; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.borrowed_tags (
    borrow_id integer NOT NULL,
    borrower_id bigint NOT NULL,
    original_tag_id integer NOT NULL,
    borrowed_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP,
    returned_at timestamp with time zone
);


ALTER TABLE public.borrowed_tags OWNER TO postgres;

--
-- Name: borrowed_tags_borrow_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.borrowed_tags_borrow_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.borrowed_tags_borrow_id_seq OWNER TO postgres;

--
-- Name: borrowed_tags_borrow_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.borrowed_tags_borrow_id_seq OWNED BY public.borrowed_tags.borrow_id;


--
-- Name: citations; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.citations (
    id integer NOT NULL,
    reference_id integer,
    user_id bigint NOT NULL,
    citation_style text NOT NULL,
    citation_text text NOT NULL,
    created_at timestamp without time zone DEFAULT now()
);


ALTER TABLE public.citations OWNER TO postgres;

--
-- Name: citations_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.citations_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.citations_id_seq OWNER TO postgres;

--
-- Name: citations_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.citations_id_seq OWNED BY public.citations.id;


--
-- Name: faction_members; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.faction_members (
    user_id bigint NOT NULL,
    faction_name character varying(255) NOT NULL
);


ALTER TABLE public.faction_members OWNER TO postgres;

--
-- Name: factions; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.factions (
    name character varying(255) NOT NULL,
    xp numeric DEFAULT 0 NOT NULL,
    level integer DEFAULT 1 NOT NULL
);


ALTER TABLE public.factions OWNER TO postgres;

--
-- Name: loop_configs; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.loop_configs (
    guild_id bigint NOT NULL,
    channel_id bigint,
    enabled boolean DEFAULT false,
    updated_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP
);


ALTER TABLE public.loop_configs OWNER TO postgres;

--
-- Name: moderation_counts; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.moderation_counts (
    user_id bigint NOT NULL,
    flagged_count integer DEFAULT 0,
    last_flagged timestamp without time zone DEFAULT now()
);


ALTER TABLE public.moderation_counts OWNER TO postgres;

--
-- Name: pdf_catalog; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.pdf_catalog (
    id integer NOT NULL,
    user_id bigint NOT NULL,
    title text NOT NULL,
    file_url text NOT NULL,
    description text,
    tags text[],
    uploaded_at timestamp without time zone DEFAULT now()
);


ALTER TABLE public.pdf_catalog OWNER TO postgres;

--
-- Name: pdf_catalog_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.pdf_catalog_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.pdf_catalog_id_seq OWNER TO postgres;

--
-- Name: pdf_catalog_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.pdf_catalog_id_seq OWNED BY public.pdf_catalog.id;


--
-- Name: pdfs; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.pdfs (
    id integer NOT NULL,
    reference_id integer,
    file_url text NOT NULL,
    created_at timestamp without time zone DEFAULT now()
);


ALTER TABLE public.pdfs OWNER TO postgres;

--
-- Name: pdfs_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.pdfs_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.pdfs_id_seq OWNER TO postgres;

--
-- Name: pdfs_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.pdfs_id_seq OWNED BY public.pdfs.id;


--
-- Name: reference_list; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.reference_list (
    id integer NOT NULL,
    user_id bigint NOT NULL,
    location_id integer NOT NULL,
    title text NOT NULL,
    authors text[] NOT NULL,
    publication_year integer,
    doi text,
    abstract text,
    tags integer[],
    created_at timestamp without time zone DEFAULT now(),
    updated_at timestamp without time zone DEFAULT now()
);


ALTER TABLE public.reference_list OWNER TO postgres;

--
-- Name: reference_list_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.reference_list_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.reference_list_id_seq OWNER TO postgres;

--
-- Name: reference_list_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.reference_list_id_seq OWNED BY public.reference_list.id;


--
-- Name: roles_backup; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.roles_backup (
    user_id bigint NOT NULL,
    role_ids bigint[],
    "timestamp" bigint
);


ALTER TABLE public.roles_backup OWNER TO postgres;

--
-- Name: tags; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.tags (
    tag_id integer NOT NULL,
    name character varying(255) NOT NULL,
    location_id bigint NOT NULL,
    owner_id bigint NOT NULL,
    content text,
    attachment_url text,
    tag_type character varying(50) DEFAULT 'default'::character varying,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP
);


ALTER TABLE public.tags OWNER TO postgres;

--
-- Name: tags_tag_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.tags_tag_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.tags_tag_id_seq OWNER TO postgres;

--
-- Name: tags_tag_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.tags_tag_id_seq OWNED BY public.tags.tag_id;


--
-- Name: users; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.users (
    id bigint NOT NULL,
    name character varying(255) NOT NULL,
    create_date timestamp with time zone DEFAULT CURRENT_TIMESTAMP,
    level integer DEFAULT 1 NOT NULL,
    exp numeric DEFAULT 0 NOT NULL,
    faction_name character varying(255) DEFAULT NULL::character varying
);


ALTER TABLE public.users OWNER TO postgres;

--
-- Name: annotations id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.annotations ALTER COLUMN id SET DEFAULT nextval('public.annotations_id_seq'::regclass);


--
-- Name: borrowed_tags borrow_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.borrowed_tags ALTER COLUMN borrow_id SET DEFAULT nextval('public.borrowed_tags_borrow_id_seq'::regclass);


--
-- Name: citations id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.citations ALTER COLUMN id SET DEFAULT nextval('public.citations_id_seq'::regclass);


--
-- Name: pdf_catalog id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdf_catalog ALTER COLUMN id SET DEFAULT nextval('public.pdf_catalog_id_seq'::regclass);


--
-- Name: pdfs id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdfs ALTER COLUMN id SET DEFAULT nextval('public.pdfs_id_seq'::regclass);


--
-- Name: reference_list id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.reference_list ALTER COLUMN id SET DEFAULT nextval('public.reference_list_id_seq'::regclass);


--
-- Name: tags tag_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.tags ALTER COLUMN tag_id SET DEFAULT nextval('public.tags_tag_id_seq'::regclass);


--
-- Data for Name: annotations; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.annotations (id, pdf_id, user_id, page_number, content, highlighted_text, created_at) FROM stdin;
\.


--
-- Data for Name: borrowed_tags; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.borrowed_tags (borrow_id, borrower_id, original_tag_id, borrowed_at, returned_at) FROM stdin;
\.


--
-- Data for Name: citations; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.citations (id, reference_id, user_id, citation_style, citation_text, created_at) FROM stdin;
\.


--
-- Data for Name: faction_members; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.faction_members (user_id, faction_name) FROM stdin;
154749533429956608	THC
\.


--
-- Data for Name: factions; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.factions (name, xp, level) FROM stdin;
THC	7.314513043214700307	1
\.


--
-- Data for Name: loop_configs; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.loop_configs (guild_id, channel_id, enabled, updated_at) FROM stdin;
\.


--
-- Data for Name: moderation_counts; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.moderation_counts (user_id, flagged_count, last_flagged) FROM stdin;
154749533429956608	3	2025-01-20 14:21:34.51824
\.


--
-- Data for Name: pdf_catalog; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.pdf_catalog (id, user_id, title, file_url, description, tags, uploaded_at) FROM stdin;
1	154749533429956608	Creatine Is a Scavenger for Methylglyoxal under Physiological Conditions via Formation of N-(4-Methyl-5-oxo-1-imidazolin-2-yl)sarcosine (MG-HCr)	/home/spawd/Downloads/pdfs/Creatine_Is_a_Scavenger_for_Methylglyoxal_under_Physiological_Conditions_via_Formation_of_N-(4-Methyl-5-oxo-1-imidazolin-2-yl)sarcosine_(MG-HCr)_154749533429956608.pdf	https://doi.org/10.1021/jf505998z	{"carbonyl stress; creatine; diabetes; dicarbonyl compounds; glycation; meat; methylglyoxal"}	2025-01-20 11:04:28.105413
\.


--
-- Data for Name: pdfs; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.pdfs (id, reference_id, file_url, created_at) FROM stdin;
\.


--
-- Data for Name: reference_list; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.reference_list (id, user_id, location_id, title, authors, publication_year, doi, abstract, tags, created_at, updated_at) FROM stdin;
\.


--
-- Data for Name: roles_backup; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.roles_backup (user_id, role_ids, "timestamp") FROM stdin;
\.


--
-- Data for Name: tags; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.tags (tag_id, name, location_id, owner_id, content, attachment_url, tag_type, created_at) FROM stdin;
\.


--
-- Data for Name: users; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.users (id, name, create_date, level, exp, faction_name) FROM stdin;
154749533429956608	spawd.	2025-03-19 17:34:03.388736-04	1	7.901556796346933191	THC
1325155980727816236	fredrick.krueger	2025-03-26 08:08:17.067277-04	1	0.695153851281538845	\N
797847699709231115	idkidkidkshauna	2025-03-26 08:09:02.02735-04	1	1.287937889035274066	\N
1286700751451848755	kartsalapsi	2025-03-26 08:09:19.128906-04	1	0.08000285962624443	\N
1036004538504708299	dinguskitty	2025-03-20 14:01:25.092254-04	1	0.048761001320651025	\N
1352971995758989373	adem08360	2025-03-24 19:27:12.768803-04	1	0.03570315333775605	\N
1353869482380234772	jubilant_dragon_20974	2025-03-24 19:28:10.130343-04	1	0.03633183023051432	\N
1353866485830783059	darine08532	2025-03-24 19:46:26.537488-04	1	0.0403606401024761	\N
832012774040141894	User_832012774040141894	2025-03-24 19:51:07.175095-04	1	0.03807772621594197	\N
1130481811965886515	beanie2722	2025-03-26 16:09:51.790644-04	1	0.159223484441386855	\N
1222656637488205866	telepathyconspiracy	2025-03-26 15:45:07.423755-04	1	0.119831615652808029	\N
\.


--
-- Name: annotations_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.annotations_id_seq', 1, false);


--
-- Name: borrowed_tags_borrow_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.borrowed_tags_borrow_id_seq', 1, false);


--
-- Name: citations_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.citations_id_seq', 1, false);


--
-- Name: pdf_catalog_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.pdf_catalog_id_seq', 1, true);


--
-- Name: pdfs_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.pdfs_id_seq', 1, false);


--
-- Name: reference_list_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.reference_list_id_seq', 1, false);


--
-- Name: tags_tag_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.tags_tag_id_seq', 1, false);


--
-- Name: annotations annotations_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.annotations
    ADD CONSTRAINT annotations_pkey PRIMARY KEY (id);


--
-- Name: borrowed_tags borrowed_tags_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.borrowed_tags
    ADD CONSTRAINT borrowed_tags_pkey PRIMARY KEY (borrow_id);


--
-- Name: citations citations_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.citations
    ADD CONSTRAINT citations_pkey PRIMARY KEY (id);


--
-- Name: faction_members faction_members_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.faction_members
    ADD CONSTRAINT faction_members_pkey PRIMARY KEY (user_id, faction_name);


--
-- Name: factions factions_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.factions
    ADD CONSTRAINT factions_pkey PRIMARY KEY (name);


--
-- Name: loop_configs loop_configs_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.loop_configs
    ADD CONSTRAINT loop_configs_pkey PRIMARY KEY (guild_id);


--
-- Name: moderation_counts moderation_counts_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.moderation_counts
    ADD CONSTRAINT moderation_counts_pkey PRIMARY KEY (user_id);


--
-- Name: pdf_catalog pdf_catalog_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdf_catalog
    ADD CONSTRAINT pdf_catalog_pkey PRIMARY KEY (id);


--
-- Name: pdfs pdfs_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdfs
    ADD CONSTRAINT pdfs_pkey PRIMARY KEY (id);


--
-- Name: reference_list reference_list_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.reference_list
    ADD CONSTRAINT reference_list_pkey PRIMARY KEY (id);


--
-- Name: roles_backup roles_backup_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.roles_backup
    ADD CONSTRAINT roles_backup_pkey PRIMARY KEY (user_id);


--
-- Name: tags tags_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.tags
    ADD CONSTRAINT tags_pkey PRIMARY KEY (tag_id);


--
-- Name: users users_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_pkey PRIMARY KEY (id);


--
-- Name: idx_borrowed_tags_unique_active; Type: INDEX; Schema: public; Owner: postgres
--

CREATE UNIQUE INDEX idx_borrowed_tags_unique_active ON public.borrowed_tags USING btree (borrower_id, original_tag_id) WHERE (returned_at IS NULL);


--
-- Name: idx_tags_lower_name_location_owner; Type: INDEX; Schema: public; Owner: postgres
--

CREATE UNIQUE INDEX idx_tags_lower_name_location_owner ON public.tags USING btree (lower((name)::text), location_id, owner_id);


--
-- Name: annotations annotations_pdf_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.annotations
    ADD CONSTRAINT annotations_pdf_id_fkey FOREIGN KEY (pdf_id) REFERENCES public.pdfs(id) ON DELETE CASCADE;


--
-- Name: borrowed_tags borrowed_tags_original_tag_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.borrowed_tags
    ADD CONSTRAINT borrowed_tags_original_tag_id_fkey FOREIGN KEY (original_tag_id) REFERENCES public.tags(tag_id) ON DELETE CASCADE;


--
-- Name: citations citations_reference_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.citations
    ADD CONSTRAINT citations_reference_id_fkey FOREIGN KEY (reference_id) REFERENCES public.reference_list(id) ON DELETE CASCADE;


--
-- Name: faction_members faction_members_faction_name_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.faction_members
    ADD CONSTRAINT faction_members_faction_name_fkey FOREIGN KEY (faction_name) REFERENCES public.factions(name) ON DELETE CASCADE;


--
-- Name: faction_members faction_members_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.faction_members
    ADD CONSTRAINT faction_members_user_id_fkey FOREIGN KEY (user_id) REFERENCES public.users(id) ON DELETE CASCADE;


--
-- Name: pdfs pdfs_reference_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdfs
    ADD CONSTRAINT pdfs_reference_id_fkey FOREIGN KEY (reference_id) REFERENCES public.reference_list(id) ON DELETE CASCADE;


--
-- Name: SCHEMA public; Type: ACL; Schema: -; Owner: pg_database_owner
--

REVOKE USAGE ON SCHEMA public FROM PUBLIC;
GRANT USAGE ON SCHEMA public TO spawd;
GRANT USAGE ON SCHEMA public TO postgres;


--
-- Name: TABLE annotations; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.annotations TO spawd;


--
-- Name: SEQUENCE annotations_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.annotations_id_seq TO spawd;


--
-- Name: TABLE borrowed_tags; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.borrowed_tags TO spawd;


--
-- Name: SEQUENCE borrowed_tags_borrow_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.borrowed_tags_borrow_id_seq TO spawd;


--
-- Name: TABLE citations; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.citations TO spawd;


--
-- Name: SEQUENCE citations_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.citations_id_seq TO spawd;


--
-- Name: TABLE faction_members; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.faction_members TO spawd;


--
-- Name: TABLE factions; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.factions TO spawd;


--
-- Name: TABLE loop_configs; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.loop_configs TO spawd;


--
-- Name: TABLE moderation_counts; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.moderation_counts TO spawd;


--
-- Name: TABLE pdf_catalog; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.pdf_catalog TO spawd;


--
-- Name: SEQUENCE pdf_catalog_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.pdf_catalog_id_seq TO spawd;


--
-- Name: TABLE pdfs; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.pdfs TO spawd;


--
-- Name: SEQUENCE pdfs_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.pdfs_id_seq TO spawd;


--
-- Name: TABLE reference_list; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.reference_list TO spawd;


--
-- Name: SEQUENCE reference_list_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.reference_list_id_seq TO spawd;


--
-- Name: TABLE roles_backup; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.roles_backup TO spawd;


--
-- Name: TABLE tags; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.tags TO spawd;


--
-- Name: SEQUENCE tags_tag_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.tags_tag_id_seq TO spawd;


--
-- Name: TABLE users; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.users TO spawd;


--
-- Name: DEFAULT PRIVILEGES FOR SEQUENCES; Type: DEFAULT ACL; Schema: public; Owner: postgres
--

ALTER DEFAULT PRIVILEGES FOR ROLE postgres IN SCHEMA public GRANT ALL ON SEQUENCES TO postgres;
ALTER DEFAULT PRIVILEGES FOR ROLE postgres IN SCHEMA public GRANT ALL ON SEQUENCES TO spawd;


--
-- Name: DEFAULT PRIVILEGES FOR FUNCTIONS; Type: DEFAULT ACL; Schema: public; Owner: postgres
--

ALTER DEFAULT PRIVILEGES FOR ROLE postgres IN SCHEMA public GRANT ALL ON FUNCTIONS TO postgres;
ALTER DEFAULT PRIVILEGES FOR ROLE postgres IN SCHEMA public GRANT ALL ON FUNCTIONS TO spawd;


--
-- Name: DEFAULT PRIVILEGES FOR TABLES; Type: DEFAULT ACL; Schema: public; Owner: postgres
--

ALTER DEFAULT PRIVILEGES FOR ROLE postgres IN SCHEMA public GRANT SELECT,INSERT,DELETE,UPDATE ON TABLES TO postgres;
ALTER DEFAULT PRIVILEGES FOR ROLE postgres IN SCHEMA public GRANT SELECT,INSERT,DELETE,UPDATE ON TABLES TO spawd;


--
-- PostgreSQL database dump complete
--

