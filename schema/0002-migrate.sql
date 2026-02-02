UPDATE active_bans
SET last_kicked = NOW()
WHERE reset = FALSE;
UPDATE active_text_mutes
SET last_muted = NOW()
WHERE reset = FALSE;
