ALTER TABLE active_caps
DROP CONSTRAINT active_caps_moderation_type_check1;
UPDATE active_caps
SET category = 'vmute'
WHERE category = 'voice_mute';

UPDATE active_caps
SET category = 'tmute'
WHERE category = 'text_mute';
ALTER TABLE active_caps
ADD CONSTRAINT active_caps_moderation_type_check
CHECK (category = ANY (ARRAY['ban','vmute','tmute']));
