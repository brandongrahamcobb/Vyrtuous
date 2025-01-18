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
    exp numeric DEFAULT 0 NOT NULL
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
-- Data for Name: loop_configs; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.loop_configs (guild_id, channel_id, enabled, updated_at) FROM stdin;
\.


--
-- Data for Name: pdf_catalog; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.pdf_catalog (id, user_id, title, file_url, description, tags, uploaded_at) FROM stdin;
4	154749533429956608	MG-HCr, the Methylglyoxal-Derived Hydroimidazolone of Creatine, a Biomarker for the Dietary Intake of Animal Source Food	https://cdn.discordapp.com/attachments/1330241956164534403/1330255483298254928/MG-HCr_the_Methylglyoxal-Derived_Hydroimidazolone_of_Creatine_a_Biomarker_for_the_Dietary_Intake_of_Animal_Source_Food_-_PubMed.pdf?ex=678d5076&is=678bfef6&hm=77a7db5847bcdd8b255f7ae8b17c17e5a147352538e569fc5251a5ae7842d13a&	https://pubs.acs.org/doi/10.1021/acs.jafc.0c00907	{methylglyoxal,glycation,creatine,"dietary study",biomarker,"meat consumption",veganism,vegetarianism}	2025-01-18 14:19:08.285982
8	154749533429956608	Plasma levels of advanced glycation end products in healthy, long-term vegetarians and subjects on a western mixed diet	https://cdn.discordapp.com/attachments/1330241956164534403/1330260537488965653/Plasma_levels_of_advanced_glycation_end_products_in_healthy_long-term_vegetarians_and_subjects_on_a_western_mixed_diet___European_Journal_of_Nutrition.pdf?ex=678d552b&is=678c03ab&hm=0c3d7c4320bff6fd78ec0ac7f9b0f47f395a23eec9fc98ca58345badc39d5526&	https://link.springer.com/article/10.1007/s394-001-8356-3	{"vegetarian diet","advanced glycation end products",carboxymethyllysine,"kindey function"}	2025-01-18 14:39:12.983434
9	154749533429956608	Morphine Binds Creatine Kinase B and Inhibits Its Activity	https://cdn.discordapp.com/attachments/1330241956164534403/1330271623689011210/fncel-12-00464.pdf?ex=678d5f7e&is=678c0dfe&hm=e8758bd11c31ea181dc2e0fd00dd83a4bd246727090b215405259bab855ffba8&	https://doi.org/10.3389/fncel.2018.00464	{ASB9,"ankyrin repeat and SOCS box protein 9; CK-B","brain creatine kinase; CK-M","muscular creatine kinase; CNS","central nervous system; I2B","I2-binding; M3G","morphine-3-glucuronide; M6G","morphine-6-glucuronide; PEBP","phosphatidylethanolamine-binding protein; SEM","standard error of the mean."}	2025-01-18 15:23:16.284516
10	154749533429956608	Effect of Creatine and Weight Training on Muscle Creatine and Performance in Vegetarians	https://cdn.discordapp.com/attachments/1330241956164534403/1330272981859045438/effect_of_creatine_and_weight_training_on_muscle.25.pdf?ex=678d60c2&is=678c0f42&hm=535189628e0df4d165c513cfb860ee4f15a99d22caea9076c9b04f44fdd3ef51&	https://journals.lww.com/acsm-msse/fulltext/2003/11000/effect_of_creatine_and_weight_training_on_muscle.25.aspx	{"LEAN TISSUE MASS","DUAL-ENERGY X-RAY ABSORPTIOMETRY","MUSCLE FIBER AREA","BIOELECTRICAL IMPEDANCE"}	2025-01-18 15:28:39.938982
11	154749533429956608	Perspective: Creatine, a Conditionally Essential Nutrient: Building the Case	https://cdn.discordapp.com/attachments/1330241956164534403/1330273450702667887/nmab111.pdf?ex=678d6131&is=678c0fb1&hm=09e4ed960a7602e284d3b2e50a9f433582d2f4945aa0faaf9dc2351e68aab181&	https://pmc.ncbi.nlm.nih.gov/articles/PMC8803499/	{creatine,food,meat,deficiency,growth,vegetarians}	2025-01-18 15:30:31.722046
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
-- Data for Name: tags; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.tags (tag_id, name, location_id, owner_id, content, attachment_url, tag_type, created_at) FROM stdin;
\.


--
-- Data for Name: users; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.users (id, name, create_date, level, exp) FROM stdin;
154749533429956608	spawd.	2025-01-18 12:34:50.695001-05	1	0.334234375096067903854191172285936772823333740234375
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

SELECT pg_catalog.setval('public.pdf_catalog_id_seq', 11, true);


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
-- Name: loop_configs loop_configs_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.loop_configs
    ADD CONSTRAINT loop_configs_pkey PRIMARY KEY (guild_id);


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
-- Name: pdfs pdfs_reference_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdfs
    ADD CONSTRAINT pdfs_reference_id_fkey FOREIGN KEY (reference_id) REFERENCES public.reference_list(id) ON DELETE CASCADE;


--
-- PostgreSQL database dump complete
--

