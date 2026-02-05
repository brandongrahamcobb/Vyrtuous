ALTER TABLE moderation_logs
ADD COLUMN executor_highest_role TEXT,
ADD COLUMN target_highest_role TEXT;
UPDATE moderation_logs
SET executor_highest_role = highest_role;
ALTER TABLE moderation_logs
DROP COLUMN highest_role;
