ALTER TABLE moderation_logs
ADD COLUMN channel_members_voice_count INTEGER DEFAULT 0,
ADD COLUMN guild_members_offline_and_online_member_count INTEGER DEFAULT 0,
ADD COLUMN guild_members_online_count INTEGER DEFAULT 0,
ADD COLUMN guild_members_voice_count INTEGER DEFAULT 0;
