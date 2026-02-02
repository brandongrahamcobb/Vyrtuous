UPDATE active_bans
SET last_kicked = NOW()
WHERE reset = FALSE;
UPDATE active_text_mutes
SET last_muted = NOW()
WHERE reset = FALSE;
DROP TRIGGER IF EXISTS set_expired_developer_logs ON bug_tracking;
