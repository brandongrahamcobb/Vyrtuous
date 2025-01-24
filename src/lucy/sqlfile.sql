--
-- PostgreSQL database dump
--

-- Dumped from database version 16.6
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

--
-- Name: update_tokens_from_pledge(); Type: FUNCTION; Schema: public; Owner: postgres
--

CREATE FUNCTION public.update_tokens_from_pledge() RETURNS trigger
    LANGUAGE plpgsql
    AS $_$
BEGIN
    UPDATE users
    SET token_balance = NEW.pledge_amount_cents / 100 * 10 -- Example: $1 = 10 tokens
    WHERE user_id = NEW.user_id;
    RETURN NEW;
END;
$_$;


ALTER FUNCTION public.update_tokens_from_pledge() OWNER TO postgres;

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
-- Name: moderation_counts; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.moderation_counts (
    user_id bigint NOT NULL,
    flagged_count integer DEFAULT 0 NOT NULL
);


ALTER TABLE public.moderation_counts OWNER TO postgres;

--
-- Name: patreon_data; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.patreon_data (
    patreon_id integer NOT NULL,
    user_id bigint,
    patreon_email character varying(255),
    pledge_amount_cents integer NOT NULL,
    tier character varying(50) NOT NULL,
    synced_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP
);


ALTER TABLE public.patreon_data OWNER TO postgres;

--
-- Name: patreon_data_patreon_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.patreon_data_patreon_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.patreon_data_patreon_id_seq OWNER TO postgres;

--
-- Name: patreon_data_patreon_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.patreon_data_patreon_id_seq OWNED BY public.patreon_data.patreon_id;


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
-- Name: reference_tags; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.reference_tags (
    reference_id bigint NOT NULL,
    tag_id bigint NOT NULL
);


ALTER TABLE public.reference_tags OWNER TO postgres;

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
-- Name: tags_tag_id_seq1; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.tags_tag_id_seq1
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.tags_tag_id_seq1 OWNER TO postgres;

--
-- Name: tags_tag_id_seq1; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.tags_tag_id_seq1 OWNED BY public.tags.tag_id;


--
-- Name: token_usage_logs; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.token_usage_logs (
    usage_id integer NOT NULL,
    user_id bigint,
    tokens_used integer NOT NULL,
    action_type character varying(255) NOT NULL,
    description text,
    created_at timestamp with time zone DEFAULT CURRENT_TIMESTAMP
);


ALTER TABLE public.token_usage_logs OWNER TO postgres;

--
-- Name: token_usage_logs_usage_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.token_usage_logs_usage_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.token_usage_logs_usage_id_seq OWNER TO postgres;

--
-- Name: token_usage_logs_usage_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.token_usage_logs_usage_id_seq OWNED BY public.token_usage_logs.usage_id;


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
-- Name: patreon_data patreon_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.patreon_data ALTER COLUMN patreon_id SET DEFAULT nextval('public.patreon_data_patreon_id_seq'::regclass);


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

ALTER TABLE ONLY public.tags ALTER COLUMN tag_id SET DEFAULT nextval('public.tags_tag_id_seq1'::regclass);


--
-- Name: token_usage_logs usage_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.token_usage_logs ALTER COLUMN usage_id SET DEFAULT nextval('public.token_usage_logs_usage_id_seq'::regclass);


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
-- Data for Name: moderation_counts; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.moderation_counts (user_id, flagged_count) FROM stdin;
1277012472016015370	1
1273029373242638366	1
1310354178882670682	1
\.


--
-- Data for Name: patreon_data; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.patreon_data (patreon_id, user_id, patreon_email, pledge_amount_cents, tier, synced_at) FROM stdin;
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
12	154749533429956608	Effects of creatine supplementation on homocysteine levels and lipid peroxidation in rats	https://cdn.discordapp.com/attachments/1330275922535780474/1330650535203700757/effects-of-creatine-supplementation-on-homocysteine-levels-and-lipid-peroxidation-in-rats.pdf?ex=678ec061&is=678d6ee1&hm=dd076dfbaa6344c72f0bad9588dfcc58cc810ebea2ac6fd1767f9582db9faeeb&	https://doi.org/10.1017/S0007114508162985	{"Creatine supplementation","Methyl balance",Homocysteine,"Lipid peroxidation"}	2025-01-19 16:30:10.012466
13	154749533429956608	Heterocyclic amines: Mutagens/carcinogens produced during cooking of meat and fish	https://cdn.discordapp.com/attachments/1330275922535780474/1330656286832070707/CAS-95-290.pdf?ex=678ec5bd&is=678d743d&hm=32cd98de7a3e534b051ab05a2a050031ed471646da6ad5f4e69df387c7545d0e&	https://doi.org/10.1111/j.1349-7006.2004.tb03205.x	\N	2025-01-19 16:53:01.480115
14	154749533429956608	Creatin(in)e and Maillard reaction products as precursors of mutagenic compounds: Effects of various amino acids	https://cdn.discordapp.com/attachments/1330275922535780474/1330658881906675753/Creatinine_and_Maillard_reaction_products_as_precursors_of_mutagenic_compounds__Effects_of_various_amino_acids___CoLab.pdf?ex=678ec827&is=678d76a7&hm=ebbce58e362f51124fece9589cbc5a161a13b12f3a6238e46b266c1aea747d6f&	https://doi.org/10.1016/0308-8146(83)90014-6	{"General Medicine","Analytical Chemistry","Food Science"}	2025-01-19 17:03:20.09157
15	154749533429956608	Creatine Is a Scavenger for Methylglyoxal under Physiological Conditions via Formation of N-(4-Methyl-5-oxo-1-imidazolin-2-yl)sarcosine (MG-HCr)	https://cdn.discordapp.com/attachments/1330275922535780474/1330663469552500756/Creatine_Is_a_Scavenger_for_Methylglyoxal_under_Physiological_Conditions_via_Formation_of_N-4-Methyl-5-oxo-1-imidazolin-2-ylsarcosine_MG-HCr___Journal_of_Agricultural_and_Food_Chemistry.pdf?ex=678ecc6d&is=678d7aed&hm=9f1394f8f0e2c3c73ea82f29b660bcc00108569adb232ceec410398efcd11f0a&	https://doi.org/10.1021/jf505998z	{"carbonyl stress; creatine; diabetes; dicarbonyl compounds; glycation; meat; methylglyoxal"}	2025-01-19 17:21:34.050811
16	154749533429956608	TESTING OF STIMULANT EFFECTS OF COFFEE ON THE PSYCHOMOTOR PERFORMANCE	https://cdn.discordapp.com/attachments/1304148357379260437/1331073500160593961/testing_of_stimulant_effects_of_coffee_on_the.3.pdf?ex=67904a4c&is=678ef8cc&hm=9fdcbf2c5ccf655d33f241af16aaa42abc8c888be65d18437bf5f21c4669de0b&	https://journals.lww.com/iphr/abstract/1997/29010/testing_of_stimulant_effects_of_coffee_on_the.3.aspx	{Coffee,"psychomotor performance","substitution test","cancellation test"}	2025-01-20 20:30:58.039917
17	154749533429956608	The Impact of Coffee on Health	https://cdn.discordapp.com/attachments/1304148357379260437/1331084439819784264/s-0043-115007.pdf?ex=6790547c&is=678f02fc&hm=5e13308196323a1851dfda9d221b421db8f08f385caf50b7850b092467fbdf4a&	https://www.thieme-connect.com/products/ejournals/pdf/10.1055/s-0043-115007.pdf	{"ALT alanine aminotransferase","CD Crohnʼs disease","CHI3L1 chitinase 3-like protein 1","CTGF connective tissue growth factor","CVD cardiovascular disease","CXCL13 small cytokine belonging the CXC family","CYP1A2 member of the cytochrome P450 oxidase system","DSS dextran sulfate sodium","GABA gamma-aminobutyric acid","GLUT 4 glucose transporter type 4","GLP 1 glucagon-like peptide 1","IBD inflammatory bowel disease","IL10 interleukin 10","NAT2 N-acetyltransferase 2","NF-κB nuclear factor kappa-light-chain-enhancer of activated B cells","NHANES National Health and Nutrition Examination Survey","TCBQ tetrachlorobenzoquinone","TNF-α tumor necrosis factor-α","TLR4 toll-like receptor 4","UC ulcerative colitis","WNT10B member of the WNT gene family"}	2025-01-20 21:14:21.077372
18	154749533429956608	The neuroprotective effects of caffeine in neurodegenerative diseases	https://cdn.discordapp.com/attachments/1304148357379260437/1331086086398541877/CNS-23-272.pdf?ex=67905605&is=678f0485&hm=a4358b9f9d4f95e91fcf021d61d9eef901fe7ffc93cdb6e1b19d1efcfbd2c15b&	https://doi.org/10.1111/cns.12684	{"adenosine receptor","Alzheimer disease","amyotrophic lateral sclerosis",caffeine,dosage,"Huntington disease","neurodegenerative disease",neuroprotection,"Parkinson disease"}	2025-01-20 21:20:53.669772
19	154749533429956608	Caffeine Synthesis and Its Mechanism and Application by Microbial Degradation, A Review	https://cdn.discordapp.com/attachments/1304148357379260437/1331090681577340928/foods-12-02721.pdf?ex=67905a4c&is=678f08cc&hm=7505528ed197058ebb0e4a323e20414271a6122ff1f8fd1ffd2b847a47258093&	https://doi.org/10.3390/foods12142721	{"caffeine; N-demethylation; C-8 oxidation; methylxanthine; microorganisms"}	2025-01-20 21:39:09.36314
20	154749533429956608	Caffeine and progression of Parkinson’s disease: A deleterious interaction with creatine	https://cdn.discordapp.com/attachments/1304148357379260437/1331094695719473152/nihms-711957.pdf?ex=67905e0a&is=678f0c8a&hm=37b057a535e88c5ccf4c858d5734dcb9e1d9973346ce9c48f698f1007d92353e&	https://doi.org/10.1097/WNF.0000000000000102	\N	2025-01-20 21:55:06.399397
21	154749533429956608	The Effect of Creatine Nitrate and Caffeine Individually or Combined on Exercise Performance and Cognitive Function: A Randomized, Crossover, Double-Blind, Placebo-Controlled Trial	https://cdn.discordapp.com/attachments/1304148357379260437/1331097490036097134/nutrients-16-00766.pdf?ex=679060a4&is=678f0f24&hm=4b007d5edd2c5a5a42d244dd15b5039afd384a6bee6c7375183fede0aaf41137&	https://doi.org/10.3390/nu16060766	{"caffeine; cognitive function; creatine nitrate; dietary supplements; ergogenic aids; exercise performance; resistance training; sports nutrition"}	2025-01-20 22:06:12.619738
22	154749533429956608	Caffeine ingestion and fluid balance: a review	https://cdn.discordapp.com/attachments/1330242076373159936/1331392742370312233/Caffeine_ingestion_and_fluid_balance__a_review_-_PubMed.pdf?ex=6791739d&is=6790221d&hm=89b00d696b09cb3fbdd02051762a125c03a801e097bdc96025ffe71c6693c98e&	https://doi.org/10.1046/j.1365-277X.2003.00477.x	\N	2025-01-21 17:39:26.322716
23	154749533429956608	Mécanismes de l’effet diurétique de la caféine	https://cdn.discordapp.com/attachments/1330242076373159936/1331394033477746799/medsci20163205p485.pdf?ex=679174d1&is=67902351&hm=758116b8fd6ab455e857ea2ee41c02cf2a9c8ca5305d730bccc7f6b741de3101&	https://doi.org/10.1051/medsci/20163205015	\N	2025-01-21 17:44:34.143634
24	154749533429956608	How empathic are vegan medical professionals compared to others? Leads from a paper–pencil-survey	https://cdn.discordapp.com/attachments/1331434796802375771/1331443954457247846/How_empathic_are_vegan_medical_professionals_compared_to_others__Leads_from_a_paperpencil-survey___European_Journal_of_Clinical_Nutrition.pdf?ex=6791a34f&is=679051cf&hm=b53d5572ef1f2f0c9fd3dceee046d3d0d9f0a8a83fb1f9fbb6c43331858ebdc3&	https://doi.org/10.1038/s41430-017-0007-8	{epidemoiology,"lifestyle modification",nutrition}	2025-01-21 21:02:56.329445
25	154749533429956608	The Combined Effect of Smoking and Coffee Drinkingon LDL and HDL Cholesterol	https://cdn.discordapp.com/attachments/1331434796802375771/1331444453503668335/heyden-et-al-1979-the-combined-effect-of-smoking-and-coffee-drinking-on-ldl-and-hdl-cholesterol.pdf?ex=6791a3c6&is=67905246&hm=e08dc53670871dec3a6410dd9a736e678f14dde89cd67ec31cceceb6398c4b5d&	https://doi.org/10.1161/01.CIR.60.1.22	\N	2025-01-21 21:04:55.21181
26	154749533429956608	The Biology of Veganism: Plasma Metabolomics Analysis Reveals Distinct Profiles of Vegans and Non-Vegetarians in the Adventist Health Study-2 Cohort	https://cdn.discordapp.com/attachments/1331434057652899940/1331445590692991017/nutrients-14-00709-v3.pdf?ex=6791a4d5&is=67905355&hm=f0938cb6e4922b7fc8ad076ab2081edebb823faf4d9f6f53b01af786e6942604&	https://doi.org/10.3390/nu14030709	\N	2025-01-21 21:09:26.375476
27	154749533429956608	Vegetarians have an indirect positive effect on sleep quality through depression condition	https://cdn.discordapp.com/attachments/1331434057652899940/1331446216243937352/s41598-023-33912-7.pdf?ex=6791a56b&is=679053eb&hm=7500616e6b329d56134c2028ffbe7d519897b5b21f4f816383c284cc26f61640&	https://doi.org/10.1038/s41598-023-33912-7	\N	2025-01-21 21:11:55.469318
28	154749533429956608	BETEL CHEWING AND CANCER.	https://cdn.discordapp.com/attachments/1331434057652899940/1331448850141151323/brmedj06288-0051c.pdf?ex=6791a7de&is=6790565e&hm=2e675f1d89e1847cb1d24638b126d1edfb460a276f0266b7ccae3bb2b13c9df8&	https://pmc.ncbi.nlm.nih.gov/articles/PMC2317373/	\N	2025-01-21 21:22:23.608388
29	154749533429956608	Psychophysiological interactions between caffeine and nicotine	https://cdn.discordapp.com/attachments/1331434796802375771/1331498386574151752/brmedj06288-0051c.pdf?ex=6791d601&is=67908481&hm=3a9e4831ef3b278bbec04895ac563a42a1d8c323ace7020bca19b98fc733dc01&	https://doi.org/10.1016/0091-3057(91)90287-C	\N	2025-01-22 00:39:13.852668
30	154749533429956608	Testosterone replacement therapy	https://cdn.discordapp.com/attachments/1304148357379260437/1331635975121801317/brmedj06288-0051c.pdf?ex=67925625&is=679104a5&hm=cd74fa7d2435ca2d75b57970dab0fb8e6cb5073ae8e0da43a65b597d065c400d&	https://doi.org/10.1111/andr.12774	{"androgen deficiency","hormonal therapy","late-onset hypogonadism",testosterone}	2025-01-22 09:45:57.568198
31	154749533429956608	S-adenosyl-L-methionine (SAM) in adults with ADHD, RS: preliminary results from an open trial	https://cdn.discordapp.com/attachments/1331463002712576074/1331644413235040286/S-adenosyl-L-methionine_SAM_in_adults_with_ADHD_RS__preliminary_results_from_an_open_trial_-_PubMed.pdf?ex=67925e00&is=67910c80&hm=718fddceb5d478ed5c78d37bf27edfdea9b97d25a46410436d1c775d2ceb14a2&	https://pubmed.ncbi.nlm.nih.gov/2236465/	\N	2025-01-22 10:19:29.541079
32	154749533429956608	Δ9-THC-Caused Synaptic and Memory Impairments Are Mediated through COX-2 Signaling	https://cdn.discordapp.com/attachments/1331939108431593527/1331939995161722910/PIIS0092867413013603.pdf?ex=67937149&is=67921fc9&hm=ff1f820e0edcae75fe4927c447870e5728b1f0490b50a68c041c7e71548a95a1&	http://doi.org/10.1016/j.cell.2013.10.042	\N	2025-01-23 05:54:01.871316
33	154749533429956608	Effects of Coffee on the Gastro-Intestinal Tract: A Narrative Review and Literature Update	https://cdn.discordapp.com/attachments/1331946271426084918/1331971310515130388/nutrients-14-00399.pdf?ex=67938e73&is=67923cf3&hm=4f346d438453296c3388f272df2a0c2d6ac43e82443187179efd01b5a7cd2a99&	https://doi.org/10.3390/nu14020399	{coffee,"gastro-intestinal tract","gastro-esophageal reflux",gallstones,"colon motility",microbiota,cancer}	2025-01-23 07:58:27.939901
34	154749533429956608	Profound effects of combining choline and piracetam on memory enhancement and cholinergic function in aged rats	https://cdn.discordapp.com/attachments/1331946271426084918/1331989754576633866/Profound_effects_of_combining_choline_and_piracetam_on_memory_enhancement_and_cholinergic_function_in_aged_rats_-_ScienceDirect.pdf?ex=67939fa0&is=67924e20&hm=0d078c3448488c61e06a56a9e15b91ed6ae0562a51a0c525ec464c335670e05d&	https://doi.org/10.1016/0197-4580(81)90007-5	\N	2025-01-23 09:11:45.351695
35	154749533429956608	Piracetam... just had the worst experience with it	https://cdn.discordapp.com/attachments/1331946271426084918/1331991263754326036/Piracetam..._just_had_the_worst_experience_with_it_-_Brain_Health_-_LONGECITY.pdf?ex=6793a108&is=67924f88&hm=7d6aeab7e5f56905d02918d4731b46f9d4f6aacc90c3cfcbed0afc3ca4f7ba1c&	https://web.archive.org/web/20230503190000/https://www.longecity.org/forum/topic/20188-piracetam-just-had-the-worst-experience-with-it/	\N	2025-01-23 09:17:45.303397
36	154749533429956608	CDP-Choline vs Alpha-GPC: Which is better for focus, memory and brain health?	https://cdn.discordapp.com/attachments/1331946271426084918/1331999945594507324/CDP-Choline_vs_Alpha-GPC__Which_is_better_for_focus_memory_and_brain_health__IJEST.pdf?ex=6793a91e&is=6792579e&hm=2b4a64adf35e0856e1d8dcf1a836b8c1d802b833f66a9efc50cce58df672e55b&	https://www.ijest.org/uncategorized/cdp-choline-vs-alpha-gpc/#	\N	2025-01-23 09:52:15.038683
37	154749533429956608	Dietary choline supplementation increases the density of nicotine binding sites in rat brain	https://cdn.discordapp.com/attachments/1331946271426084918/1332009892789227520/document.pdf?ex=6793b261&is=679260e1&hm=0fce14e099f7dec3d974a2d6dae864c010dcce3c0d93b46f83662871964f1ced&	https://pubmed.ncbi.nlm.nih.gov/1326624/	\N	2025-01-23 10:31:46.622295
38	154749533429956608	Effects of Prenatal Nicotine Exposure on Primate Brain Development and Attempted Amelioration with Supplemental Choline or Vitamin C: Neurotransmitter Receptors, Cell Signaling and Cell Development Biomarkers in Fetal Brain Regions of Rhesus Monkeys	https://cdn.discordapp.com/attachments/1331946271426084918/1332014330547277845/Effects_of_Prenatal_Nicotine_Exposure_on_Primate_B.pdf?ex=6793b684&is=67926504&hm=ffa1a5e969383ae7224616642ece92bfc5d01b786573e5505f00a7f51f8b2100&	http://doi.org/10.1038/sj.npp.1300544	\N	2025-01-23 10:49:24.754488
39	154749533429956608	Caffeine potentiates the enhancement by choline of striatal acetylcholine release	https://cdn.discordapp.com/attachments/1331946271426084918/1332014585196187688/Caffeine_potentiates_the_enhancement_by_choline_of_striatal_acetylcholine_release_-_PubMed.pdf?ex=6793b6c0&is=67926540&hm=dce93ef0708a00d565332404db3e2ba2b9088e7b179d03e58ad7da910359f322&	https://doi.org/10.1016/0024-3205(92)90622-v	\N	2025-01-23 10:50:25.371268
40	154749533429956608	Phytochemical overview and medicinal importance of Coffea species from the past until now	https://cdn.discordapp.com/attachments/1331946271426084918/1332021908635516968/1-s2.0-S1995764516304680-main.pdf?ex=6793bd92&is=67926c12&hm=5dffb55fa9648778c02e3833e35d1f020b54a6aa8e2e94c3b7ba2266de685276&	https://doi.org/10.1016/j.apjtm.2016.11.008	{Phytochemistry,"Medicinal importance",Caffeine,Polyphenols}	2025-01-23 11:19:31.520052
41	154749533429956608	Effect of Supercritical Carbon Dioxide Decaffeination on Volatile Components of Green Teas	https://cdn.discordapp.com/attachments/1331946271426084918/1332023900133130353/Effect_of_Supercritical_Carbon_Dioxide_Decaffeination_on_Volatile_Components_of_Green_Teas_-_Lee_-_2007_-_Journal_of_Food_Science_-_Wiley_Online_Library.pdf?ex=6793bf6d&is=67926ded&hm=28a072bda151f73fc2ca807356cbdab75cb484d2fbac4a4aa08e5ca76d16c97c&	https://doi.org/10.1111/j.1750-3841.2007.00446.x	\N	2025-01-23 11:27:26.460937
42	154749533429956608	Total Water-Soluble Choline Concentration Does Not Differ in Milk from Vegan, Vegetarian, and Nonvegetarian Lactating Women	https://cdn.discordapp.com/attachments/1331946271426084918/1332026737785311305/1-s2.0-S0022316622020594-main.pdf?ex=6793c212&is=67927092&hm=b761fd689c3e1b2bedd51fd9f8c27867d664c96e65aef0f04dec4838375e377b&	https://doi.org/10.1093/jn/nxz257	{choline,"breast milk",vegan,vegetarian,phosphocholine,glycerophosphocholine}	2025-01-23 11:38:42.826253
43	154749533429956608	Vitamin B-12 content in breast milk of vegan, vegetarian, and nonvegetarian lactating women in the United States	https://cdn.discordapp.com/attachments/1331946271426084918/1332028759196438579/1-s2.0-S0002916522029471-mainext.pdf?ex=6793c3f4&is=67927274&hm=c2d07410eee63c08df5bb19c4ce4b56ad1cb52096775bc98ea42eaafa3e00ce2&	https://doi.org/10.1093/ajcn/nqy104	{"vitamin B-12",vegan,vegetarian,"breast milk",supplements}	2025-01-23 11:46:44.925233
44	154749533429956608	Dietary Reference Values for choline	https://cdn.discordapp.com/attachments/1331946271426084918/1332034174315663492/EFSA_Journal_-_2016_-_-_Dietary_Reference_Values_for_choline.pdf?ex=6793c8ff&is=6792777f&hm=e6a373a1751ff1824f3d9a738edba53315530f43ecf99967581fd60b40c4e553&	https://doi.org/10.2903/j.efsa.2016.4484	{choline,phosphatidylcholine,"observed intake","depletion/repletion study","Adequate Intake","Dietary Reference Value"}	2025-01-23 12:08:15.998471
45	154749533429956608	The Relationship between Choline Bioavailability from Diet, Intestinal Microbiota Composition, and Its Modulation of Human Diseases	https://cdn.discordapp.com/attachments/1331946071005335615/1332037132604149901/nutrients-12-02340.pdf?ex=6793cbc0&is=67927a40&hm=a8d43634af9e4de30ea18236025d7db6650e776364f13a8a72030eb8e0bfb872&	https://doi.org/10.3390/nu12082340	{choline,TMA,TMAO,"non-alcoholic steatohepatitis (NASH)","cardiovascular disease (CVD)","chronic kidney diseases (CKD)",probiotics,"gut microbiota",polyphenols,"fecal microbiota transplantation"}	2025-01-23 12:20:01.121865
46	154749533429956608	The combined effects of L-theanine and caffeine on cognitive performance and mood	https://cdn.discordapp.com/attachments/1331946071005335615/1332078870874034238/The_combined_effects_of_L-theanine_and_caffeine_on_cognitive_performance_and_moo.pdf?ex=6793f29f&is=6792a11f&hm=edaebe58a9f4c9f779ff966184379841ec82a7c54edfed98fa4b42878542332e&	https://doi.org/10.1179/147683008X301513	\N	2025-01-23 15:05:52.376012
47	154749533429956608	The neuropharmacology of L-theanine(N-ethyl-L-glutamine): a possible neuroprotective and cognitive enhancing agent	https://cdn.discordapp.com/attachments/1331946071005335615/1332122436103180390/The_neuropharmacology_of_L-theanineN-ethyl-L-glutamine-_a_possible_neuroprotec.pdf?ex=67941b32&is=6792c9b2&hm=c82ca04c08da7e8d22429adbce3de72460f1e2dcc0e56bc236d108a6bb352954&	https://pubmed.ncbi.nlm.nih.gov/17182482/	\N	2025-01-23 17:58:59.142439
48	154749533429956608	L-Theanine Prevents Long-Term Affective and Cognitive Side Effects of Adolescent D-9-Tetrahydrocannabinol\nExposure and Blocks Associated Molecular and Neuronal Abnormalities in the Mesocorticolimbic Circuitry	https://cdn.discordapp.com/attachments/1331946071005335615/1332123352034447500/zns739.pdf?ex=67941c0c&is=6792ca8c&hm=1fef1d66e6d70aa84d7a048ecad126df01d00bfa91c8b3a31143d7de66d52912&	https://doi.org/10.1523/JNEUROSCI.1050-20.2020	{adolescence,cognition,dopamine,l-theanine,mesocorticolimbic,THC}	2025-01-23 18:02:37.557997
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
-- Data for Name: reference_tags; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.reference_tags (reference_id, tag_id) FROM stdin;
\.


--
-- Data for Name: tags; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.tags (tag_id, name, location_id, owner_id, content, attachment_url, tag_type, created_at) FROM stdin;
\.


--
-- Data for Name: token_usage_logs; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.token_usage_logs (usage_id, user_id, tokens_used, action_type, description, created_at) FROM stdin;
\.


--
-- Data for Name: users; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.users (id, name, create_date, level, exp) FROM stdin;
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
-- Name: patreon_data_patreon_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.patreon_data_patreon_id_seq', 1, false);


--
-- Name: pdf_catalog_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.pdf_catalog_id_seq', 48, true);


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
-- Name: tags_tag_id_seq1; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.tags_tag_id_seq1', 1, false);


--
-- Name: token_usage_logs_usage_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.token_usage_logs_usage_id_seq', 1, false);


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
-- Name: moderation_counts moderation_counts_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.moderation_counts
    ADD CONSTRAINT moderation_counts_pkey PRIMARY KEY (user_id);


--
-- Name: patreon_data patreon_data_patreon_email_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.patreon_data
    ADD CONSTRAINT patreon_data_patreon_email_key UNIQUE (patreon_email);


--
-- Name: patreon_data patreon_data_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.patreon_data
    ADD CONSTRAINT patreon_data_pkey PRIMARY KEY (patreon_id);


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
-- Name: reference_tags reference_tags_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.reference_tags
    ADD CONSTRAINT reference_tags_pkey PRIMARY KEY (reference_id, tag_id);


--
-- Name: tags tags_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.tags
    ADD CONSTRAINT tags_pkey PRIMARY KEY (tag_id);


--
-- Name: token_usage_logs token_usage_logs_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.token_usage_logs
    ADD CONSTRAINT token_usage_logs_pkey PRIMARY KEY (usage_id);


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
-- Name: patreon_data trigger_update_tokens; Type: TRIGGER; Schema: public; Owner: postgres
--

CREATE TRIGGER trigger_update_tokens AFTER INSERT OR UPDATE ON public.patreon_data FOR EACH ROW EXECUTE FUNCTION public.update_tokens_from_pledge();


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
-- Name: TABLE loop_configs; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.loop_configs TO spawd;


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
-- Name: TABLE tags; Type: ACL; Schema: public; Owner: postgres
--

GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE public.tags TO spawd;


--
-- Name: SEQUENCE tags_tag_id_seq; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.tags_tag_id_seq TO spawd;


--
-- Name: SEQUENCE tags_tag_id_seq1; Type: ACL; Schema: public; Owner: postgres
--

GRANT ALL ON SEQUENCE public.tags_tag_id_seq1 TO spawd;


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

