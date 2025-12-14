INSERT INTO statistic_channels (
    channel_id,
    enabled,
    guild_id,
    snowflakes,
    type
)
SELECT
    channel_id,
    COALESCE(enabled, FALSE),
    guild_id,
    snowflakes,
    COALESCE(type, 'general')
FROM log_channels
ON CONFLICT (guild_id, channel_id) DO NOTHING;

ALTER TABLE users
RENAME COLUMN server_muter_guild_ids TO administrator_guild_ids
DELETE COLUMN coordinator_room_names;
DELETE COLUMN moderator_room_names;
