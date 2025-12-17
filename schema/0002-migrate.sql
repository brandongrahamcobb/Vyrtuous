ALTER TABLE active_caps
ADD COLUMN duration_seconds INTEGER;
ALTER TABLE active_caps
DROP COLUMN duration;

