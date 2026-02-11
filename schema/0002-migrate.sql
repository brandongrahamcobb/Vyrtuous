ALTER TABLE public.command_aliases
DROP CONSTRAINT command_aliases_category_check;

ALTER TABLE public.command_aliases
ADD CONSTRAINT command_aliases_category_check
CHECK (category = ANY (ARRAY['vegan','vmute','ban','flag','tmute','role']));
DROP TABLE ban_roles;
DROP TABLE hide_roles;
DROP TABLE users;
DROP TABLE text_mute_roles;
DROP SEQUENCE IF EXISTS public.history_id_seq CASCADE;

DROP FUNCTION IF EXISTS public.set_expired() CASCADE;
DROP TRIGGER IF EXISTS set_expired_active_bans ON public.active_bans;
DROP TRIGGER IF EXISTS set_expired_active_stages ON public.active_stages;
DROP TRIGGER IF EXISTS set_expired_active_text_mutes ON public.active_text_mutes;
DROP TRIGGER IF EXISTS set_expired_active_voice_mutes ON public.active_voice_mutes;
