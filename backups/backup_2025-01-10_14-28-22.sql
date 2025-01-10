--
-- PostgreSQL database dump
--

-- Dumped from database version 16.6
-- Dumped by pg_dump version 16.6

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
-- Name: loop_configs; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.loop_configs (
    guild_id bigint NOT NULL,
    channel_id bigint NOT NULL,
    enabled boolean DEFAULT false NOT NULL
);


ALTER TABLE public.loop_configs OWNER TO postgres;

--
-- Name: tags; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.tags (
    id integer NOT NULL,
    name character varying(100) NOT NULL,
    location_id bigint NOT NULL,
    content text,
    attachment_url text,
    owner_id bigint NOT NULL,
    tag_type character varying(10) NOT NULL,
    CONSTRAINT tags_tag_type_check CHECK (((tag_type)::text = ANY (ARRAY[('default'::character varying)::text, ('loop'::character varying)::text])))
);


ALTER TABLE public.tags OWNER TO postgres;

--
-- Name: tags_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.tags_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.tags_id_seq OWNER TO postgres;

--
-- Name: tags_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.tags_id_seq OWNED BY public.tags.id;


--
-- Name: tags id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.tags ALTER COLUMN id SET DEFAULT nextval('public.tags_id_seq'::regclass);


--
-- Data for Name: loop_configs; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.loop_configs (guild_id, channel_id, enabled) FROM stdin;
730907954345279591	985926652041261117	t
1300517536001036348	1326623116536844308	t
\.


--
-- Data for Name: tags; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.tags (id, name, location_id, content, attachment_url, owner_id, tag_type) FROM stdin;
1	bookClub	730907954345279591	📚 **Book Recommendation:** ***'The Joyful Vegan: How to Stay Vegan in a World That Wants You to Eat Meat, Dairy, and Egg'* by Colleen Patrick-Goudreau**\n\n*In these pages, Colleen shares her wisdom for managing these challenges and arms readers—both vegan and plant-based—with solutions and strategies for 'coming out vegan' to family, friends, and colleagues; cultivating healthy relationships (with vegans and non-vegans); communicating effectively; sharing enthusiasm without proselytizing; finding like-minded community; and experiencing peace of mind as a vegan in a non-vegan world.  \n\nBy implementing the tools provided in this book, readers will find they can live ethically, eat healthfully, engage socially—and remain a joyful vegan.*	\N	154749533429956608	loop
2	discordOutreach	730907954345279591	Want a chance to do activism through **Discord Outreach** with us in ARA? Discord Outreach is an activism event we host weekly in the Events tab where activists will be assigned into groups to join the VCs of larger servers to do vegan activism. \n🔹An Outreach Leader will lead the group by coordinating which server and which VC channel the group will join.\n🔹Groups will have speakers and listeners. Speakers will mainly direct activism discussion for nonvegan users who join the channel. Listeners will take up spots in the channel to prevent other random users from joining to prevent distraction from ongoing discussion.\n🔹You can sign up to be pinged for Discord Outreach in <#990761562199457813>	\N	154749533429956608	loop
3	tipsVegan	730907954345279591	**Looking for tips to meet local vegans or activism events? **\n🔹Try searching for vegan Facebook groups for your closet major city or area. \n🔹Get in contact with an animal rights organization like [PETA](<https://www.peta.org/>), [Direct Action Everywhere](<https://www.directactioneverywhere.com/>), [Mercy for Animals](<https://mercyforanimals.org/>), [Humane Society of the US](<https://www.humanesociety.org/>), etc in your area. Try searching for an organization promoting plant-based eating as well! \n🔹You can also search in [Meetup](<https://www.meetup.com/home/>), a social media platform for organizing events and activities. \n🔹Volunteering at animal sanctuaries.\n🔹Start a Facebook or [Meetup](<https://www.meetup.com/home/>) group yourself!	\N	154749533429956608	loop
4	restaurantVegan	730907954345279591	**Been to a nice vegan or vegan-friendly restaurant recently? **\n🔹Write a Google or Yelp review with the keyword 'vegan' about your experience and what foods you had.\n🔹Share your opinion on the restaurant on Instagram or Facebook, bonus if you post in vegan/vegetarian Facebook groups.\n🔹Sign up for the **[Happy Cow](<https://www.happycow.net/>)** 🐮 💜  app to list the restaurant and leave a review to help others easily access veg options too.	\N	154749533429956608	loop
5	happyCow	730907954345279591	Sign up for the **[Happy Cow](<https://www.happycow.net/>)** 🐮 💜 app, a mobile app and website that lists vegan and vegan-friendly restaurants and also a passionate community of over one million vegan-focused members. Aside from listing restaurants it also lists farmers markets, health food stores and all types of businesses with a vegan focus.	\N	154749533429956608	loop
6	veganProduct	730907954345279591	**Enjoyed a vegan product recently? **\n🔹Share your opinion on the product on Instagram or Facebook, bonus if you post in vegan Facebook groups.\n🔹Sign up on the **[abillion](<https://www.abillion.com/>)** app and write your review of the product. The platform allows users to review plant-based, cruelty-free and sustainable products, while donating between 0.10 and $1 to nonprofit organizations for each review written.	\N	154749533429956608	loop
7	vegRecipe	730907954345279591	**Tried out a great online recipe recently?** Be sure to leave a high rating and review to boost your favorite vegan and plant-based creator! ⭐	\N	154749533429956608	loop
8	localOutreach	730907954345279591	**Looking for street outreach opportunities?** Try searching for any local chapters from [Anonymous for the Voiceless](<https://www.anonymousforthevoiceless.org/>), [We The Free](<https://www.activism.wtf/>), or events in vegan Facebook/[Meetup](<https://www.meetup.com/home/>) groups.	\N	154749533429956608	loop
9	getPolitical	730907954345279591	**Get political!** Join in local pressure campaigns and getting ballot measures passed with groups such as [Animal Activist Mentorship](<https://www.animalactivismmentorship.com/>), [PETA](<https://www.peta.org/action/campaigns/>), [Plant Based Treaty](<https://plantbasedtreaty.org/>), & [Pro-Animal Future](<https://proanimal.org/>) in the US, [Viva!](<https://viva.org.uk/>) in the UK, [Animal Justice Party](<https://www.animaljusticeparty.org/>) in AU. They can be coalitions with goals ranging from banning fur, banning foi gras, banning cages, getting plant-based milks in schools, to banning factory farms.	\N	154749533429956608	loop
10	nonConActivism	730907954345279591	**Prefer a more non-confrontational form of activism?** Consider these ideas!\n🔹Sidewalk chalking is a great way for public visual messaging in public areas with higher foot traffic. Chalking is not permanent and non-damaging so legally does not typically count as vandalism, so it's usually allowed but check with your local municipals first.\n🔹You can also draw vegan messaging at the beach in the sand, weather permitting.\n🔹Consider using vegan usernames like in online gaming or social media.\n🔹Stickering such as placing them on your phone, car, laptop, or water bottle when going out.\n🔹Wearing clothes that promote the vegan message whether you're just going out for groceries or at the gym to show off your cruelty-free gains to others.	\N	154749533429956608	loop
11	veganFacebook	730907954345279591	**Be a foundation for local vegan community building.**\n🔹There are numerous Facebook groups to assist new vegans and the veg curious in finding resources in their local community. If you don't have one, consider starting one yourself!\n🔹Direct people to **[r/Vegan](<https://www.reddit.com/r/vegan/>)** or **[r/AskVegans](<https://www.reddit.com/r/AskVegans/>)** on Reddit to ask questions or utilize the search function in the groups for specific advice.\n🔹Schedule vegan potlucks, game nights, or other events on [Meetup](<https://www.meetup.com/home/>) for your area.	\N	154749533429956608	loop
12	veganSkills	730907954345279591	**Harness your skills!** Utilize your unique skills and talent to be in service for the animals, such as:\n🔹If you are a programmer or software engineer, consider volunteering with [Vegan Hacktivists](<https://veganhacktivists.org/>).\n🔹If you are a graphic designer, you can help design pamphlets or T-shirts.\n🔹If you are handy, consider volunteering at animal sanctuaries to help construct infrastructure for the residents.\n🔹If you're a cook, consider taking photos and posting them in social media or foodie groups. \n🔹If you got music or comedic talent, consider going to open mic events or volunteering at Veg Fests about veganism.	\N	154749533429956608	loop
13	betterOutreach	730907954345279591	**How can I become a better Outreacher?**\n🔹Here is a useful **[video](<https://www.youtube.com/watch?v=-nznQXhXgMY>)** on a conversation structure guide by [The Victim's Perspective](<https://www.youtube.com/@TheVictimsPerspective>) on Youtube\n🔹Learn from prominent vegan outreachers like [Earthling Ed](<https://www.youtube.com/@ed.winters/featured>), [Joey Carbstrong](<https://www.youtube.com/@JoeyCarbstrong>), [Debug Your Brain](<https://www.youtube.com/@DebugYourBrain>), [Clif Grant](<https://www.youtube.com/@clifgrant>), [David Ramms](<https://www.youtube.com/@davidramms>), and more by watching their content. \n🔹Do group outreach with [Anonymous for the Voiceless](<https://www.anonymousforthevoiceless.org/>), [We The Free](<https://www.activism.wtf/>), or host activism events in vegan Facebook/[Meetup](<https://www.meetup.com/home/>) groups.	\N	154749533429956608	loop
14	onlineComment	730907954345279591	**Online comment section activism ideas:**\n🔹Leaving comments on viral videos or posts on veganism or related videos that can direct toward veganism.\n🔹Getting vegan allies involved to help give a 'like' to your comment or post to get noticed. \n🔹Carnists giving you a short fuse? Consider keeping a digital document with saved pre-written replies to copy and paste to help avoid being tempted to use condescending tone in replies.	\N	154749533429956608	loop
15	getActive	730907954345279591	**Willing to your get your hands dirty and be proactive for the animals?**\n🔹Consider doing direct action, attending vigils such as by the [Animal Save Movement](<https://thesavemovement.org/>), or get involved in pressure campaigns.\n🔹Get in contact with organizations such as [Direct Action Everywhere](<https://www.directactioneverywhere.com/>), [Animal Rebellion](<https://animalrebellion.org/about/>), [PETA](<https://www.peta.org/action/campaigns/>), and [Animal Liberation Front](<https://animalliberationfrontline.com/>)\n🔹**[Video](<https://www.youtube.com/watch?v=LHyqJxSeUFc>)** on the importance of pressure campaigns by [The Cranky Vegan](<https://www.instagram.com/the.cranky.vegan/?hl=en>) on [VeganFTA](<https://veganfta.com/>)	\N	154749533429956608	loop
16	USAVegan	730907954345279591	**In the US?** Get into legislation activism! Find state and local representatives to send letters about animal rights, meat subsidies, ag gag laws, environmental impacts, or increased food disease risks to by using [CommonCause.org](<https://www.commoncause.org/find-your-representative/>) by entering your street address.	\N	154749533429956608	loop
17	USA2Vegan	730907954345279591	**In the US 🇺🇸?** Get connected with [Agriculture Fairness Alliance](<https://www.agriculturefairnessalliance.org/>) for legislation activism! A 501(c)(4) nonprofit whose mission is to strategically employ lobbyists to accelerate policy changes that make sustainable plant-based food accessible to everyone at a price they can afford, empower communities to develop local plant based agriculture systems, and give farmers tools and strategy to transition from animal ag to plant based farming.	\N	154749533429956608	loop
18	DEVegan	730907954345279591	**In Germany 🇩🇪?** Get connected with [V-Party3](<https://v-partei.de/>) for legislation activism!\n*Die V-Partei ist eine deutsche Partei, die der Tierproduktindustrie den Kampf angesagt hat, mit Verboten jeglicher tierischen Produkten, Tierversuchen und Zurschaustellung in Zoo und Circus. Zusätzlich setzen sie sich für ernstzunehmende ethische und ernährungstechnische Bildung an Schulen und den Schutz von Tierrechtsaktivisten ein und fördern bezahlbare Nahrungsmittel aus solidarischer Landwirtschaft mit ökologischen Alternativen zu Pestiziden.*	\N	154749533429956608	default
19	cheatsheet	730907954345279591	**Need a cheatsheet for responding to justifications to harm and exploit animals?** Bookmark [Vegan Sidekick](<https://www.godfist.com/vegansidekick/guide.php>). A comprehensive list of the known excuses and the responses for them. It will also link to the common comebacks after responding to certain excuses too!	\N	154749533429956608	loop
20	stream	730907954345279591	**Looking for a streaming service of Plant-Based News & Entertainment Network for FREE? **Download [UnchainedTV](<https://unchainedtv.com/>) on your phone via the APP store or on your TV via your Amazon Fire Stick, AppleTV device or Roku device.	\N	154749533429956608	loop
21	contentCreator	730907954345279591	**Need some food content creator recommendations? **\n🔹[Nora Cooks](<https://www.noracooks.com/>) - Recipes that are easy to make and even easier to eat\n🔹[Forks Over Knives](<https://www.forksoverknives.com/recipes/>) - Healthy whole food plant-based recipes\n🔹[Rainbow Plant Life](<https://rainbowplantlife.com/>) - For the home cook looking to wow their friends\n🔹[Cheap Lazy Vegan](<https://thecheaplazyvegan.com/blog/>) - Easy and affordable vegan meal ideas\n🔹[The Foodie Takes Flight](<https://thefoodietakesflight.com/>) - Asian-inspired recipes\n🔹[Vegan Richa](<https://www.veganricha.com/recipes/>) - Indian-inspired recipes\n🔹[Eat Figs Not Pigs](<https://www.eatfigsnotpigs.com/>) - Fusion comfort foods\n🔹[Thee Burger Dude](<https://theeburgerdude.com/>) - Popular fast food recipes veganized	\N	154749533429956608	loop
22	nutrition	730907954345279591	**Need a comprehensive source on vegan nutrition?** [Vegan Health](<https://veganhealth.org/>) is a website with sources and studies by registered dieticians on evidence-based nutrient recommendations.	\N	154749533429956608	loop
23	calendar	730907954345279591	**Looking for activism opportunities and events near you?** Try [Animal Rights Calender](<https://animalrightscalendar.com/>)! Not finding an organization or event? Contact email: *person@animalrightscalendar.com*	\N	154749533429956608	loop
24	gallery	730907954345279591	**Need a database of stock animal rights images for activism?** Try <https://stock.weanimals.org/>	\N	154749533429956608	loop
27	horseBook	730907954345279591	📚 **Book Recommendation:** ***'Riding On the Power of Others: A Horsewoman's Path to Unconditional Love'* by Ren Hurst**\n\n*Ren Hurst's memoir explores her journey of self-discovery and healing through her relationships with horses, learning to let go of control and embrace unconditional love and acceptance that leads her to walk away completely from riding and training horses and into a world where relationship is all that matters.. \n\nThrough her experiences, Hurst reveals the transformative power of horses and the natural world to teach humans about empathy, compassion, and the interconnectedness of all beings.*	\N	154749533429956608	loop
28	restaurantVegan	730907954345279591	\N	\N	259318132936671233	default
\.


--
-- Name: tags_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.tags_id_seq', 28, true);


--
-- Name: loop_configs loop_configs_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.loop_configs
    ADD CONSTRAINT loop_configs_pkey PRIMARY KEY (guild_id);


--
-- Name: tags tags_name_location_id_owner_id_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.tags
    ADD CONSTRAINT tags_name_location_id_owner_id_key UNIQUE (name, location_id, owner_id);


--
-- Name: tags tags_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.tags
    ADD CONSTRAINT tags_pkey PRIMARY KEY (id);


--
-- PostgreSQL database dump complete
--

