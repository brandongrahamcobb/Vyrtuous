--
-- PostgreSQL database dump
--

-- Dumped from database version 17.4
-- Dumped by pg_dump version 17.4

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
154749533429956608	Dream
\.


--
-- Data for Name: factions; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.factions (name, xp, level) FROM stdin;
Dream	0.81217959977011187	1
\.


--
-- Data for Name: loop_configs; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.loop_configs (guild_id, channel_id, enabled, updated_at) FROM stdin;
1347284827350630591	1347284828894265398	t	2025-04-05 14:31:46.876127-04
\.


--
-- Data for Name: moderation_counts; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.moderation_counts (user_id, flagged_count, last_flagged) FROM stdin;
858565178979385354	7	2025-04-01 07:27:29.420417
791945428832223242	12	2025-04-04 12:49:10.918224
154749533429956608	0	2025-01-20 14:21:34.51824
366138954497523725	2	2025-04-04 18:56:20.945411
\.


--
-- Data for Name: pdf_catalog; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.pdf_catalog (id, user_id, title, file_url, description, tags, uploaded_at) FROM stdin;
1	154749533429956608	Creatine Is a Scavenger for Methylglyoxal under Physiological Conditions via Formation of N-(4-Methyl-5-oxo-1-imidazolin-2-yl)sarcosine (MG-HCr)	/home/spawd/Downloads/pdfs/Creatine_Is_a_Scavenger_for_Methylglyoxal_under_Physiological_Conditions_via_Formation_of_N-(4-Methyl-5-oxo-1-imidazolin-2-yl)sarcosine_(MG-HCr)_154749533429956608.pdf	https://doi.org/10.1021/jf505998z	{"carbonyl stress; creatine; diabetes; dicarbonyl compounds; glycation; meat; methylglyoxal"}	2025-01-20 11:04:28.105413
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
49	154749533429956608	Effect of organic acids on the solid-state polymorphic phase transformation of piracetam	https://cdn.discordapp.com/attachments/1332517637456138343/1333244901672091689/Effect_of_organic_acids_on_the_solid-state_polymorphic_phase_transformation_of_piracetam_-_ScienceDirect.pdf?ex=67983093&is=6796df13&hm=92a1a202ea534b121ca12e413b71200798dd6ceb29ec8ba6e8a1ead9676f509a&	https://doi.org/10.1016/j.ijpharm.2023.123532	{"Crystallization kinetics; Metastable polymorph; Organic acids; Piracetam; Solid-state phase transformation."}	2025-01-26 20:19:16.064697
50	154749533429956608	Stability Hierarchy between Piracetam Forms I, II, and III from Experimental Pressure-Temperature Diagrams and Topological Inferences	https://cdn.discordapp.com/attachments/1332517637456138343/1333246652252684359/Stability_Hierarchy_between_Piracetam_Forms.pdf?ex=67983234&is=6796e0b4&hm=15403bf3acd516820d759e9b8164f64b040d2193c32902a19b7c5604e484e2e4&	https://doi.org/10.1016/j.ijpharm.2015.11.036	{"Crystal polymorphism","phase diagram","phase transition","physical stability",preformulation,"solid state"}	2025-01-26 20:26:13.499111
51	154749533429956608	Piracetam Defines a New Binding Site for Allosteric Modulators of α-Amino-3-hydroxy-5-methyl-4-isoxazole-propionic Acid (AMPA) Receptors	https://cdn.discordapp.com/attachments/1332517637456138343/1333251415019356180/Piracetam_Defines_a_New_Binding_Site_for_Allosteric_Modulators_of_-Amino-3-hydroxy-5-methyl-4-isoxazole-propionic_Acid_AMPA_Receptors___Journal_of_Medicinal_Chemistry.pdf?ex=679836a3&is=6796e523&hm=d957d529c4f03a2c530168fa6d68c93b873356beb87098d2ddb0a7d93fa5dd70&	https://doi.org/10.1021/jm901905j	{"Chemical structure",modulators,oligomers,receptors,"screening assays"}	2025-01-26 20:45:09.008728
52	154749533429956608	Improved mitochondrial function in brain aging and Alzheimer disease – the new mechanism of action of the old metabolic enhancer piracetam	https://cdn.discordapp.com/attachments/1332517637456138343/1333253776181825546/fnins-04-00044.pdf?ex=679838d6&is=6796e756&hm=201429fd2cba71117387ea516f6852fe10e17032c1a4dc89041606634f37d914&	https://doi.org/10.3389/fnins.2010.00044	{"mitochondrial dysfunction","alzheimer’s disease",aging,"oxidative stress",piracetam}	2025-01-26 20:54:31.908395
53	154749533429956608	Clinical efficacy of piracetam in cognitive impairment: a meta-analysis	https://cdn.discordapp.com/attachments/1332517637456138343/1333255476770508930/Clinical_Efficacy_of_Piracetam_in_Cognitive_Impairment__A_Meta-Analysis___Dementia_and_Geriatric_Cognitive_Disorders___Karger_Publishers.pdf?ex=67983a6c&is=6796e8ec&hm=f9acb7756de3d0817d629cb1fdcc9c995365199c7855413205077983e7a5ac17&	https://doi.org/10.1159/000057700	{Piracetam,"cognition disorders",meta-analysis,"randomised clinical trials"}	2025-01-26 21:01:17.41746
54	154749533429956608	Piracetam improves mitochondrial dysfunction following oxidative stress	https://cdn.discordapp.com/attachments/1332517637456138343/1333257315377348668/147-0706459a.pdf?ex=67983c22&is=6796eaa2&hm=a9031d5ccdc542234a2f02c26e392f776921d9e1d0ed1107e12f909c6f1a0bde&	https://doi.org/10.1038/sj.bjp.0706459	{"Piracetam; mitochondrial function; oxidative stress; aging",AD,"Alzheimer’s disease; AMPA","a-amino-3-hydroxy-5-methyl-4-isoxazole-propionic acid; ANOVA","analysis of variance; CNS","central nervous system; GPx","glutathione peroxidase; GR","glutathione reductase; H2O2","hydrogen peroxide; NaN3","sodium azide; PC12 cells","rat pheochromocytoma cells; R123","rhodamine 123; ROS","reactive oxygen species; SNP","sodium nitroprusside; SOD","superoxide dismutase; TMRE","tetramethylrhodamineethylester; TTFA",thenoyltrifluoroacetone}	2025-01-26 21:08:35.696859
55	154749533429956608	The appetite-suppressant effect of nicotine is enhanced by caffeine	https://cdn.discordapp.com/attachments/1332517637456138343/1333370027025629194/The_appetitesuppressant_effect_of_nicotine_is_enhanced_by_caffeine__-_Jessen_-_2005_-_Diabetes_Obesity_and_Metabolism_-_Wiley_Online_Library.pdf?ex=6798a51b&is=6797539b&hm=bae0e6d967e69f7e99f9e31fd4a9830cf691f7ebf8b9bd88fab58f61d5d03c92&	https://doi.org/10.1111/j.1463-1326.2004.00389.x	\N	2025-01-27 04:36:28.344112
56	154749533429956608	The Role of β-alanine Supplementation on Muscle Carnosine and Exercise Performance	https://cdn.discordapp.com/attachments/1332517637456138343/1333380339082723359/The_Role_of__alanine_Supplementation_on.pdf?ex=6798aeb5&is=67975d35&hm=4fbde9216d03d6cb45e676f39b729d85aa54f25b7a6077db415bd6c23db69a2d&	https://d1wqtxts1xzle7.cloudfront.net/43844955/The_Role_of_-alanine_Supplementation_on_20160317-20037-pvj2h0-libre.pdf?1458282741=&response-content-disposition=inline%3B+filename%3DThe_Role_of__alanine_Supplementation_on.pdf	{"buffer capacity; fatigue; physical capacity; skeletal muscle"}	2025-01-27 05:17:26.877723
57	154749533429956608	International society of sports nutrition position stand: Beta-Alanine	https://cdn.discordapp.com/attachments/1332517637456138343/1333487322297274399/s12970-015-0090-y.pdf?ex=67991258&is=6797c0d8&hm=18fa4f38b0c1ba7cf00fce5fafe5d366a369836da436518bd61c12bf49dd1e8d&	https://doi.org/10.1186/s12970-015-0090-y	{Carnosine,"Exercise Bout","Ergogenic Effect",Muscle,Carnosine,"Neuromuscular Fatigue"}	2025-01-27 12:22:33.798781
58	154749533429956608	The basics of phosphate metabolism	https://cdn.discordapp.com/attachments/1332517637456138343/1333808003106738267/gfad188.pdf?ex=679a3d00&is=6798eb80&hm=c5450570a9dd1aaf9cb288affd1448986a74170bf5958b8de2124e1fde60d772&	https://doi.org/10.1093/ndt/gfad188	{bone,"cell metabolism","endocrine regulation",intestine,kidney}	2025-01-28 09:36:50.146038
59	154749533429956608	Role of creatine and phosphocreatine in neuronal protection from anoxic and ischemic damage	https://cdn.discordapp.com/attachments/1332517637456138343/1333943313581211750/Role_of_creatine_and_phosphocreatine_in_neuronal_protection_from_anoxic_and_ischemic_damage___Amino_Acids.pdf?ex=679abb05&is=67996985&hm=5f2bdc4746937f4454a9a6385b7004604df5795fecfba9d11a63b9d614ba4d14&	https://doi.org/10.1007/s00726-001-0133-3	\N	2025-01-28 18:34:30.605649
60	154749533429956608	Creatine in the brain	https://cdn.discordapp.com/attachments/1332517637456138343/1333944498258051083/6_215.pdf?ex=679abc1f&is=67996a9f&hm=005c3a800efba53101da9d2e12301dd084c3539ad35a2b77509f37fda39612cc&	https://www.jstage.jst.go.jp/article/jpfsm/6/4/6_215/_pdf	{creatine,cognition,"magnetic resonance spectroscopy"}	2025-01-28 18:39:13.021102
61	154749533429956608	The Assertive Brain: Anterior Cingulate Phosphocreatine plus Creatine Levels Correlate With Self-Directedness in Healthy Adolescents	https://cdn.discordapp.com/attachments/1332517637456138343/1333945833611202585/fpsyt-10-00763.pdf?ex=679abd5e&is=67996bde&hm=08845bf327a6a2643d293c3a50aa5f676322097f753962e9f918efe657ae1b10&	https://www.frontiersin.org/journals/psychiatry/articles/10.3389/fpsyt.2019.00763/full	{"magnetic resonance spectroscopy","temperament character inventory",adolescence,"brain biochemistry","brain metabolism"}	2025-01-28 18:44:31.340277
62	154749533429956608	Effects of the transdermal nicotine patch on normalization of HDL-C and its subfractions	https://cdn.discordapp.com/attachments/1332517637456138343/1334328810493181962/Effects_of_the_Transdermal_Nicotine_Patch_on_Normalization_of_HDL-C_and_Its_Subfractions_-_ScienceDirect.pdf?ex=679c220b&is=679ad08b&hm=4b638dbb01d6052a4dd8d8385a69374868e0b90a4b8a485ec52d17b45dc5da3e&	https://doi.org/10.1006/pmed.2000.0694	\N	2025-01-29 20:06:20.283863
63	154749533429956608	Neuronal Ca2+ Channels and Nicotinic ACh Receptors as Functional Targets of the Nootropic Nefiracetam	https://cdn.discordapp.com/attachments/1332517637456138343/1334329753213206568/Psychogeriatrics_-_2007_-_Yoshii_-_Neuronal_Ca2_Channels_and_Nicotinic_ACh_Receptors_as_Functional_Targets_of_the.pdf?ex=679c22eb&is=679ad16b&hm=2f29e36ee724e40a2379971a21b87593c63c3f9c87ea83bf9c5af83a5105b2ef&	https://doi.org/10.1111/j.1479-8301.2001.tb00071.x	{nefeiracetam,nootropics,"neuronal Ca2+ channels","presynaptic nicotinic ACh receptor","Alzheimer’s disease"}	2025-01-29 20:10:04.958741
64	154749533429956608	Nicotine: From Discovery to Biological Effects	https://cdn.discordapp.com/attachments/1332517637456138343/1334330270605639711/ijms-24-14570.pdf?ex=679c2367&is=679ad1e7&hm=413454372d06b6fdfd48163f63c75e352209a895c864d2ef5a0840590a5d94cc&	https://doi.org/10.3390/ijms241914570	\N	2025-01-29 20:12:08.34921
65	154749533429956608	Neuronal Nicotinic Acetylcholine Receptor Structure and Function and Response to Nicotine	https://cdn.discordapp.com/attachments/1332517637456138343/1334550144594739261/nihms766536.pdf?ex=679cf02d&is=679b9ead&hm=854b1aa0ca6535a3ca048d3f26ffd4a2ba13879ce29734fe71a58875e895b316&	https://doi.org/10.1016/bs.irn.2015.07.001	\N	2025-01-30 10:45:50.59433
66	154749533429956608	Nicotinic Acetylcholine Receptors and Nicotine Addiction: A Brief Introduction	https://cdn.discordapp.com/attachments/1332517637456138343/1334564807407439912/nihms-1618249.pdf?ex=679cfdd5&is=679bac55&hm=88bca601b9b97ce816f23f24979d3f8af1c61c97172069fb6b7e7a5d776c5420&	https://doi.org/10.1016/j.neuropharm.2020.108256	\N	2025-01-30 11:44:06.451062
67	154749533429956608	Sweet Taste Signaling: The Core Pathways and Regulatory Mechanisms	https://cdn.discordapp.com/attachments/1332517637456138343/1336723845088018524/ijms-23-08225.pdf?ex=67a4d897&is=67a38717&hm=972c6bd00bd8cb918901496eb4273c73a4941a1ba9e5bd6ff7e3ffe454b9a56d&	https://doi.org/10.3390/ijms23158225	{"sweet taste receptor; gustation; G-protein-coupled receptor; noncaloric sweeteners"}	2025-02-05 10:43:21.598931
68	154749533429956608	Acetylcholine as a neuromodulator: cholinergic signaling shapes nervous system function and behavior\n	https://cdn.discordapp.com/attachments/1336333713159491614/1339046349303447632/nihms405384.pdf?ex=67ad4b97&is=67abfa17&hm=94ec09d6135749cdd2d70347c153b9796d5a54d002535ea13499201ddf5da588&	https://pmc.ncbi.nlm.nih.gov/articles/PMC3466476/	\N	2025-02-11 20:32:10.228315
69	154749533429956608	Preventing Illegitimate Extrasynaptic Acetylcholine Receptor Clustering Requires the RSU-1 Protein	https://cdn.discordapp.com/attachments/1332517637456138343/1339273356779978862/6525.full.pdf?ex=67ae1f02&is=67accd82&hm=81b28efea603a6498b8b8fc2b91ae49c4dc0b6830833d4f1b68b42e37b24931f&	https://doi.org/10.1523/JNEUROSCI.3733-15.2016	{"acetylcholine receptor","C. elegans","forward genetic screen","neuromuscular junction",RSU-1,synapse}	2025-02-12 11:34:13.026317
70	154749533429956608	Role of perisynaptic parameters in neurotransmitter homeostasis - computational study of a general synapse	https://cdn.discordapp.com/attachments/1332517637456138343/1339301870002245704/nihms666759.pdf?ex=67ae3990&is=67ace810&hm=bf17a6d629e1c91ede76e02b15b2671288d6a8b1a15e46dbc1c7456449c1b5ac&	https://doi.org/10.1002/syn.21547	{"computational model","neurotransmitter homeostasis","glial configurations","non-synaptic sources",synapse}	2025-02-12 13:27:31.062571
71	154749533429956608	Effects of soy isoflavones on cognitive function: a systematic review and meta-analysis of randomized controlled trials	https://cdn.discordapp.com/attachments/1332517637456138343/1340823804989014137/nuz050.pdf?ex=67b3c2fa&is=67b2717a&hm=8c57e7a08b004ee93f8f361db6339245f15eb0cb95f28e54da9b7aa88ce4f0a0&	https://doi.org/10.1093/nutrit/nuz050	{cognition,isoflavones,meta-analysis}	2025-02-16 18:15:09.05978
72	154749533429956608	The multiple roles of the α7 nicotinic acetylcholine receptor in modulating glutamatergic systems in the normal and diseased nervous system	https://cdn.discordapp.com/attachments/1332517637456138343/1341537802591932516/The_multiple_roles_of_the_7_nicotinic_acetylcholine_receptor_in_modulating_glutamatergic_systems_in_the_normal_and_diseased_nervous_system_-_ScienceDirect.pdf?ex=67b65bf0&is=67b50a70&hm=5595fa57d2677a59cd4a5a55ff027f2b38a204e001f3ed958fdf050278a89ee2&	https://doi.org/10.1016/j.bcp.2015.07.018	\N	2025-02-18 17:32:19.399151
73	154749533429956608	Acute Subjective Effects of Psychedelics within and Beyond WEIRD Contexts	https://cdn.discordapp.com/attachments/1332517637456138343/1343360127226351646/Acute_Subjective_Effects_of_Psychedelics_within_and_Beyond_WEIRD_Contexts_-_PubMed.pdf?ex=67bcfd1c&is=67bbab9c&hm=6ea2c22cb9028bccd9d17078b40a27e06f6dfb010288fef4cfb730328319e151&	https://doi.org/10.1080/02791072.2023.2255274	{"Psychedelics; acute psychedelic subjective effects; cross-cultural; mystical-type experience"}	2025-02-23 18:13:35.6625
74	154749533429956608	Astrocyte regulation of synaptic signaling in psychiatric disorders	https://cdn.discordapp.com/attachments/1332517637456138343/1344681326237716531/s41386-022-01338-w.pdf?ex=67c1cb93&is=67c07a13&hm=d043a2c1d175c2f00baf07e315f55eaf029e345b7f9961376542528356f9cb74&	https://doi.org/10.1038/s41386-022-01338-w	\N	2025-02-27 09:43:34.355725
75	154749533429956608	Solubility Rules for Ionic Compounds	https://cdn.discordapp.com/attachments/1343669857949847652/1345816817121235105/Solubility_Rules_for_Ionic_Compounds.pdf?ex=67c5ed15&is=67c49b95&hm=eb3daa9164ca89ad6c6b77cb45100ea954e5933071a14d47685a5f96cc3b14aa&	https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/materials-science-and-engineering/solid-state-synthesis/solubility-rules-solubility-of-common-ionic-compounds	\N	2025-03-02 12:55:36.825132
76	154749533429956608	Inhibition of Recombinant Human T-type Calcium Channels by Δ9-Tetrahydrocannabinol and Cannabidiol	https://cdn.discordapp.com/attachments/1347284828894265398/1348002276941369344/16124.pdf?ex=67cde073&is=67cc8ef3&hm=74da7500e67ce71fe7c3a41cd3eb6ac3041a8f803809aa3374e1b36b924e7e9e&	https://pmc.ncbi.nlm.nih.gov/articles/PMC3259625/	\N	2025-03-08 13:39:51.882537
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
1	bookClub	730907954345279591	730907954345279591	📚 **Book Recommendation:** ***'The Joyful Vegan: How to Stay Vegan in a World That Wants You to Eat Meat, Dairy, and Egg'* by Colleen Patrick-Goudreau**\n\n*In these pages, Colleen shares her wisdom for managing these challenges and arms readers—both vegan and plant-based—with solutions and strategies for 'coming out vegan' to family, friends, and colleagues; cultivating healthy relationships (with vegans and non-vegans); communicating effectively; sharing enthusiasm without proselytizing; finding like-minded community; and experiencing peace of mind as a vegan in a non-vegan world.  \n\nBy implementing the tools provided in this book, readers will find they can live ethically, eat healthfully, engage socially—and remain a joyful vegan.*	bookClub	default	2025-04-06 14:15:48.000703-04
2	discordOutreach	730907954345279591	730907954345279591	Want a chance to do activism through **Discord Outreach** with us in ARA? Discord Outreach is an activism event we host weekly in the Events tab where activists will be assigned into groups to join the VCs of larger servers to do vegan activism. \n🔹An Outreach Leader will lead the group by coordinating which server and which VC channel the group will join.\n🔹Groups will have speakers and listeners. Speakers will mainly direct activism discussion for nonvegan users who join the channel. Listeners will take up spots in the channel to prevent other random users from joining to prevent distraction from ongoing discussion.\n🔹You can sign up to be pinged for Discord Outreach in <#990761562199457813>	discordOutreach	default	2025-04-06 14:15:48.000703-04
3	tipsVegan	730907954345279591	730907954345279591	**Looking for tips to meet local vegans or activism events? **\n🔹Try searching for vegan Facebook groups for your closet major city or area. \n🔹Get in contact with an animal rights organization like [PETA](<https://www.peta.org/>), [Direct Action Everywhere](<https://www.directactioneverywhere.com/>), [Mercy for Animals](<https://mercyforanimals.org/>), [Humane Society of the US](<https://www.humanesociety.org/>), etc in your area. Try searching for an organization promoting plant-based eating as well! \n🔹You can also search in [Meetup](<https://www.meetup.com/home/>), a social media platform for organizing events and activities. \n🔹Volunteering at animal sanctuaries.\n🔹Start a Facebook or [Meetup](<https://www.meetup.com/home/>) group yourself!	tipsVegan	default	2025-04-06 14:15:48.000703-04
5	happyCow	730907954345279591	730907954345279591	Sign up for the **[Happy Cow](<https://www.happycow.net/>)** 🐮 💜 app, a mobile app and website that lists vegan and vegan-friendly restaurants and also a passionate community of over one million vegan-focused members. Aside from listing restaurants it also lists farmers markets, health food stores and all types of businesses with a vegan focus.	happyCow	default	2025-04-06 14:15:48.000703-04
6	veganProduct	730907954345279591	730907954345279591	**Enjoyed a vegan product recently? **\n🔹Share your opinion on the product on Instagram or Facebook, bonus if you post in vegan Facebook groups.\n🔹Sign up on the **[abillion](<https://www.abillion.com/>)** app and write your review of the product. The platform allows users to review plant-based, cruelty-free and sustainable products, while donating between 0.10 and $1 to nonprofit organizations for each review written.	veganProduct	default	2025-04-06 14:15:48.000703-04
7	vegRecipe	730907954345279591	730907954345279591	**Tried out a great online recipe recently?** Be sure to leave a high rating and review to boost your favorite vegan and plant-based creator! ⭐	vegRecipe	default	2025-04-06 14:15:48.000703-04
8	localOutreach	730907954345279591	730907954345279591	**Looking for street outreach opportunities?** Try searching for any local chapters from [Anonymous for the Voiceless](<https://www.anonymousforthevoiceless.org/>), [We The Free](<https://www.activism.wtf/>), or events in vegan Facebook/[Meetup](<https://www.meetup.com/home/>) groups.	localOutreach	default	2025-04-06 14:15:48.000703-04
9	getPolitical	730907954345279591	730907954345279591	**Get political!** Join in local pressure campaigns and getting ballot measures passed with groups such as [Animal Activist Mentorship](<https://www.animalactivismmentorship.com/>), [PETA](<https://www.peta.org/action/campaigns/>), [Plant Based Treaty](<https://plantbasedtreaty.org/>), & [Pro-Animal Future](<https://proanimal.org/>) in the US, [Viva!](<https://viva.org.uk/>) in the UK, [Animal Justice Party](<https://www.animaljusticeparty.org/>) in AU. They can be coalitions with goals ranging from banning fur, banning foi gras, banning cages, getting plant-based milks in schools, to banning factory farms.	getPolitical	default	2025-04-06 14:15:48.000703-04
10	nonConActivism	730907954345279591	730907954345279591	**Prefer a more non-confrontational form of activism?** Consider these ideas!\n🔹Sidewalk chalking is a great way for public visual messaging in public areas with higher foot traffic. Chalking is not permanent and non-damaging so legally does not typically count as vandalism, so it's usually allowed but check with your local municipals first.\n🔹You can also draw vegan messaging at the beach in the sand, weather permitting.\n🔹Consider using vegan usernames like in online gaming or social media.\n🔹Stickering such as placing them on your phone, car, laptop, or water bottle when going out.\n🔹Wearing clothes that promote the vegan message whether you're just going out for groceries or at the gym to show off your cruelty-free gains to others.	nonConActivism	default	2025-04-06 14:15:48.000703-04
11	veganFacebook	730907954345279591	730907954345279591	**Be a foundation for local vegan community building.**\n🔹There are numerous Facebook groups to assist new vegans and the veg curious in finding resources in their local community. If you don't have one, consider starting one yourself!\n🔹Direct people to **[r/Vegan](<https://www.reddit.com/r/vegan/>)** or **[r/AskVegans](<https://www.reddit.com/r/AskVegans/>)** on Reddit to ask questions or utilize the search function in the groups for specific advice.\n🔹Schedule vegan potlucks, game nights, or other events on [Meetup](<https://www.meetup.com/home/>) for your area.	veganFacebook	default	2025-04-06 14:15:48.000703-04
12	veganSkills	730907954345279591	730907954345279591	**Harness your skills!** Utilize your unique skills and talent to be in service for the animals, such as:\n🔹If you are a programmer or software engineer, consider volunteering with [Vegan Hacktivists](<https://veganhacktivists.org/>).\n🔹If you are a graphic designer, you can help design pamphlets or T-shirts.\n🔹If you are handy, consider volunteering at animal sanctuaries to help construct infrastructure for the residents.\n🔹If you're a cook, consider taking photos and posting them in social media or foodie groups. \n🔹If you got music or comedic talent, consider going to open mic events or volunteering at Veg Fests about veganism.	veganSkills	default	2025-04-06 14:15:48.000703-04
13	betterOutreach	730907954345279591	730907954345279591	**How can I become a better Outreacher?**\n🔹Here is a useful **[video](<https://www.youtube.com/watch?v=-nznQXhXgMY>)** on a conversation structure guide by [The Victim's Perspective](<https://www.youtube.com/@TheVictimsPerspective>) on Youtube\n🔹Learn from prominent vegan outreachers like [Earthling Ed](<https://www.youtube.com/@ed.winters/featured>), [Joey Carbstrong](<https://www.youtube.com/@JoeyCarbstrong>), [Debug Your Brain](<https://www.youtube.com/@DebugYourBrain>), [Clif Grant](<https://www.youtube.com/@clifgrant>), [David Ramms](<https://www.youtube.com/@davidramms>), and more by watching their content. \n🔹Do group outreach with [Anonymous for the Voiceless](<https://www.anonymousforthevoiceless.org/>), [We The Free](<https://www.activism.wtf/>), or host activism events in vegan Facebook/[Meetup](<https://www.meetup.com/home/>) groups.	betterOutreach	default	2025-04-06 14:15:48.000703-04
14	onlineComment	730907954345279591	730907954345279591	**Online comment section activism ideas:**\n🔹Leaving comments on viral videos or posts on veganism or related videos that can direct toward veganism.\n🔹Getting vegan allies involved to help give a 'like' to your comment or post to get noticed. \n🔹Carnists giving you a short fuse? Consider keeping a digital document with saved pre-written replies to copy and paste to help avoid being tempted to use condescending tone in replies.	onlineComment	default	2025-04-06 14:15:48.000703-04
15	getActive	730907954345279591	730907954345279591	**Willing to your get your hands dirty and be proactive for the animals?**\n🔹Consider doing direct action, attending vigils such as by the [Animal Save Movement](<https://thesavemovement.org/>), or get involved in pressure campaigns.\n🔹Get in contact with organizations such as [Direct Action Everywhere](<https://www.directactioneverywhere.com/>), [Animal Rebellion](<https://animalrebellion.org/about/>), [PETA](<https://www.peta.org/action/campaigns/>), and [Animal Liberation Front](<https://animalliberationfrontline.com/>)\n🔹**[Video](<https://www.youtube.com/watch?v=LHyqJxSeUFc>)** on the importance of pressure campaigns by [The Cranky Vegan](<https://www.instagram.com/the.cranky.vegan/?hl=en>) on [VeganFTA](<https://veganfta.com/>)	getActive	default	2025-04-06 14:15:48.000703-04
16	USAVegan	730907954345279591	730907954345279591	**In the US?** Get into legislation activism! Find state and local representatives to send letters about animal rights, meat subsidies, ag gag laws, environmental impacts, or increased food disease risks to by using [CommonCause.org](<https://www.commoncause.org/find-your-representative/>) by entering your street address.	USAVegan	default	2025-04-06 14:15:48.000703-04
17	USA2Vegan	730907954345279591	730907954345279591	**In the US 🇺🇸?** Get connected with [Agriculture Fairness Alliance](<https://www.agriculturefairnessalliance.org/>) for legislation activism! A 501(c)(4) nonprofit whose mission is to strategically employ lobbyists to accelerate policy changes that make sustainable plant-based food accessible to everyone at a price they can afford, empower communities to develop local plant based agriculture systems, and give farmers tools and strategy to transition from animal ag to plant based farming.	USA2Vegan	default	2025-04-06 14:15:48.000703-04
18	DEVegan	730907954345279591	730907954345279591	**In Germany 🇩🇪?** Get connected with [V-Party3](<https://v-partei.de/>) for legislation activism!\n*Die V-Partei ist eine deutsche Partei, die der Tierproduktindustrie den Kampf angesagt hat, mit Verboten jeglicher tierischen Produkten, Tierversuchen und Zurschaustellung in Zoo und Circus. Zusätzlich setzen sie sich für ernstzunehmende ethische und ernährungstechnische Bildung an Schulen und den Schutz von Tierrechtsaktivisten ein und fördern bezahlbare Nahrungsmittel aus solidarischer Landwirtschaft mit ökologischen Alternativen zu Pestiziden.*	DEVegan	default	2025-04-06 14:15:48.000703-04
19	cheatsheet	730907954345279591	730907954345279591	**Need a cheatsheet for responding to justifications to harm and exploit animals?** Bookmark [Vegan Sidekick](<https://www.godfist.com/vegansidekick/guide.php>). A comprehensive list of the known excuses and the responses for them. It will also link to the common comebacks after responding to certain excuses too!	cheatsheet	default	2025-04-06 14:15:48.000703-04
20	stream	730907954345279591	730907954345279591	**Looking for a streaming service of Plant-Based News & Entertainment Network for FREE? **Download [UnchainedTV](<https://unchainedtv.com/>) on your phone via the APP store or on your TV via your Amazon Fire Stick, AppleTV device or Roku device.	stream	default	2025-04-06 14:15:48.000703-04
21	contentCreator	730907954345279591	730907954345279591	**Need some food content creator recommendations? **\n🔹[Nora Cooks](<https://www.noracooks.com/>) - Recipes that are easy to make and even easier to eat\n🔹[Forks Over Knives](<https://www.forksoverknives.com/recipes/>) - Healthy whole food plant-based recipes\n🔹[Rainbow Plant Life](<https://rainbowplantlife.com/>) - For the home cook looking to wow their friends\n🔹[Cheap Lazy Vegan](<https://thecheaplazyvegan.com/blog/>) - Easy and affordable vegan meal ideas\n🔹[The Foodie Takes Flight](<https://thefoodietakesflight.com/>) - Asian-inspired recipes\n🔹[Vegan Richa](<https://www.veganricha.com/recipes/>) - Indian-inspired recipes\n🔹[Eat Figs Not Pigs](<https://www.eatfigsnotpigs.com/>) - Fusion comfort foods\n🔹[Thee Burger Dude](<https://theeburgerdude.com/>) - Popular fast food recipes veganized	contentCreator	default	2025-04-06 14:15:48.000703-04
22	nutrition	730907954345279591	730907954345279591	**Need a comprehensive source on vegan nutrition?** [Vegan Health](<https://veganhealth.org/>) is a website with sources and studies by registered dieticians on evidence-based nutrient recommendations.	nutrition	default	2025-04-06 14:15:48.000703-04
23	calendar	730907954345279591	730907954345279591	**Looking for activism opportunities and events near you?** Try [Animal Rights Calender](<https://animalrightscalendar.com/>)! Not finding an organization or event? Contact email: *person@animalrightscalendar.com*	calendar	default	2025-04-06 14:15:48.000703-04
24	gallery	730907954345279591	730907954345279591	**Need a database of stock animal rights images for activism?** Try <https://stock.weanimals.org/>	gallery	default	2025-04-06 14:15:48.000703-04
27	horseBook	730907954345279591	730907954345279591	📚 **Book Recommendation:** ***'Riding On the Power of Others: A Horsewoman's Path to Unconditional Love'* by Ren Hurst**\n\n*Ren Hurst's memoir explores her journey of self-discovery and healing through her relationships with horses, learning to let go of control and embrace unconditional love and acceptance that leads her to walk away completely from riding and training horses and into a world where relationship is all that matters.. \n\nThrough her experiences, Hurst reveals the transformative power of horses and the natural world to teach humans about empathy, compassion, and the interconnectedness of all beings.*	horseBook	default	2025-04-06 14:15:48.000703-04
28	restaurantVegan	730907954345279591	730907954345279591	\N	restaurantVegan	default	2025-04-06 14:15:48.000703-04
29	restaurantVegan	730907954345279591	730907954345279591	\N	restaurantVegan	default	2025-04-06 14:15:48.000703-04
31	DiscordOutreach	730907954345279591	730907954345279591	Want a chance to do activism through **Discord Outreach** with us in ARA? Discord Outreach is an activism event we host weekly in the Events tab where activists will be assigned into groups to join the VCs of larger servers to do vegan activism. \n🔹An Outreach Leader will lead the group by coordinating which server and which VC channel the group will join.\n🔹Groups will have speakers and listeners. Speakers will mainly direct activism discussion for nonvegan users who join the channel. Listeners will take up spots in the channel to prevent other random users from joining to prevent distraction from ongoing discussion.\n🔹You can sign up to be pinged for Discord Outreach in <#990761562199457813>	DiscordOutreach	default	2025-04-06 14:15:48.000703-04
1	src	1347284827350630591	154749533429956608	Vyrtuous Source Code	https://cdn.discordapp.com/attachments/1347284828894265398/1358507379490296100/lucy-16.8.0-py3-none-any.whl?ex=67f41817&is=67f2c697&hm=9391e2fe9136668a6207d52bf0b20e73df9a48df78e813ab611d126a06d92684&	default	2025-04-06 14:19:04.212568-04
2	1P-LSD	1347284827350630591	154749533429956608	\N	https://cdn.discordapp.com/attachments/1347284828894265398/1358739918087651339/IMG_5648.png?ex=67f4f0a8&is=67f39f28&hm=2f760cd10ce8a4da2900b9d46e4000295d867feb37262e846ef8f57e081176e4&	default	2025-04-07 05:43:04.81308-04
\.


--
-- Data for Name: users; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.users (id, name, create_date, level, exp, faction_name) FROM stdin;
791945428832223242	kenerion	2025-04-04 12:49:10.200468-04	1	0.704428415057684379	\N
797847699709231115	idkidkidkshauna	2025-03-26 08:09:02.02735-04	1	1.287937889035274066	\N
1286700751451848755	kartsalapsi	2025-03-26 08:09:19.128906-04	1	0.129028514396365565	\N
858565178979385354	mastermedic1929	2025-04-01 06:52:10.397861-04	1	6.060064093083376721	\N
1352971995758989373	adem08360	2025-03-24 19:27:12.768803-04	1	0.03570315333775605	\N
1353869482380234772	jubilant_dragon_20974	2025-03-24 19:28:10.130343-04	1	0.03633183023051432	\N
1353866485830783059	darine08532	2025-03-24 19:46:26.537488-04	1	0.0403606401024761	\N
832012774040141894	User_832012774040141894	2025-03-24 19:51:07.175095-04	1	0.03807772621594197	\N
1130481811965886515	beanie2722	2025-03-26 16:09:51.790644-04	1	0.159223484441386855	\N
1222656637488205866	telepathyconspiracy	2025-03-26 15:45:07.423755-04	1	0.119831615652808029	\N
1325155980727816236	fredrick.krueger	2025-03-26 08:08:17.067277-04	1	1.123347501373537534	\N
977391770734301255	myaboo_	2025-03-30 20:45:49.100722-04	1	0.087682678357341305	\N
154749533429956608	spawd.	2025-03-19 17:34:03.388736-04	3	86.308266933333634850	Dream
1036004538504708299	dinguskitty	2025-03-20 14:01:25.092254-04	1	0.124259990351212255	\N
366138954497523725	vampii	2025-04-04 15:05:33.188256-04	1	3.212666701434011291	\N
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

SELECT pg_catalog.setval('public.tags_tag_id_seq', 2, true);


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
-- Name: users users_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_pkey PRIMARY KEY (id);


--
-- Name: idx_borrowed_tags_unique_active; Type: INDEX; Schema: public; Owner: postgres
--

CREATE UNIQUE INDEX idx_borrowed_tags_unique_active ON public.borrowed_tags USING btree (borrower_id, original_tag_id) WHERE (returned_at IS NULL);


--
-- Name: annotations annotations_pdf_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.annotations
    ADD CONSTRAINT annotations_pdf_id_fkey FOREIGN KEY (pdf_id) REFERENCES public.pdfs(id) ON DELETE CASCADE;


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

