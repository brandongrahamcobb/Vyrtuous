
from vyrtuous.service.check_service import *
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.member_service import MemberService
from vyrtuous.utils.duration import Duration
from vyrtuous.utils.statistics import Statistics
from vyrtuous.utils.vegans import Vegans
from vyrtuous.utils.emojis import Emojis

import datetime

class Aliases(commands.Cog):

    def __init__(self, bot: DiscordBot):
        self.bot = bot
        self.emoji = Emojis()
        self.channel_service = ChannelService()
        self.member_service = ChannelService()
        self.statistics = Statistics()
        self.vegans = Vegans.get_vegans()

    async def handle_ban_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None
        duration = args[1] if len(args) > 1 else '24h'
        reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.'
        #'Tag a member or include their snowflake ID'),
        #'(+|-)duration(m|h|d) \n 0 - permanent / 24h - default \n `+` to append, `-` to delete, `=` to overwrite reason'),
        #'Optional reason (required for 7 days or more)')
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        allowed = False
        updated_reason = reason
        now = datetime.now(timezone.utc)
        updated_reason = reason
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj:
            if member_obj.id in self.vegans:
                return await message.reply(content=f'\U0001F6AB You cannot ban a superhero.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot ban the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB No permission to use `{alias_name}` in {channel_obj.mention}.')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB Cannot ban `{highest_role}` in {channel_obj.mention}.')
        is_coordinator = await is_owner_developer_coordinator_via_alias(message, alias)
        async with self.bot.db_pool.acquire() as conn:
            existing_ban = await conn.fetchrow('''
                SELECT expires_at, reason
                FROM active_bans
                WHERE guild_id = $1
                    AND discord_snowflake = $2
                    AND channel_id = $3
                    AND room_name = $4
            ''', message.guild.id, member_obj.id, channel_obj.id, channel_obj.name)
            duration_obj = Duration()
            reason_obj = Reason()
            duration_obj.load_from_combined_duration_str(duration)
            if existing_ban and duration not in ['+', '=', '-']:
                expires_at = existing_ban['expires_at']
                old_reason = existing_ban['reason']
                reason_obj.load_old_reason(old_reason)
            if duration not in ['+', '=', '-']:
                expires_at = duration_obj.output_datetime()
                action = None
            else:
                reason_obj.load_new_reason(reason)
                action = reason_obj.interpret_action(duration_obj.get_prefix())
                if action in ('delete', 'overwrite') and not is_coordinator:
                    return await message.reply(content='\U0001F6AB Only coordinators are allowed to overwrite the reason.')
                updated_reason = reason_obj.output_display()
            if expires_at and expires_at <= now:
                return await message.reply(content='\U0001F6AB You cannot reduce a ban below the current time.')
            duration_obj.load_base(expires_at)
            duration_display = duration_obj.output_display()
            caps = await self.get_caps_for_channel(message.guild.id, channel_obj.id)
            active_cap = next((c for c in caps if c[0] == 'ban'), None)
            if active_cap:
                cap_expires_at = active_cap[1]
            else:
                cap_expires_at = timedelta(days=7) + now
            if existing_ban and not expires_at < existing_ban['expires_at'] and not (is_coordinator or expires_at <= cap_expires_at):
                return await message.reply(content='\U0001F6AB Only coordinators can ban for longer than the channel cap.')
        try:
            await channel_obj.set_permissions(member_obj, view_channel=False, reason=f'{updated_reason}')
        except discord.Forbidden:
            return await message.reply(content=f'\U0001F6AB {member_obj.mention} was not successfully banned.', allowed_mentions=discord.AllowedMentions.none())
        is_in_channel = False
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            is_in_channel = True
            try:
                await member_obj.move_to(None, reason=f'{updated_reason}')
            except discord.Forbidden:
                await message.reply(content=f'\U0001F6AB Could not disconnect {member_obj.mention} from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            except Exception as e:
                logger.exception(f'Unexpected error while disconnecting user: {e}')
                raise
        async with self.bot.db_pool.acquire() as conn:
            await conn.execute('''
                INSERT INTO active_bans (guild_id, discord_snowflake, channel_id, reason, expires_at, room_name)
                VALUES ($1,$2,$3,$4,$5,$6)
                ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name) DO UPDATE
                SET reason=$4, expires_at=$5
            ''', message.guild.id, member_obj.id, channel_obj.id, updated_reason, expires_at, channel_obj.name)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1, $2, $3, $4, $5, $6)
            ''', 'ban', member_obj.id, message.author.id, message.guild.id, channel_obj.id, updated_reason)
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} has been banned",
            description=f"**User:** {member_obj.mention}\n**Channel:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason}",
            color=discord.Color.orange()
        )
        await message.reply(embed=embed, allowed_mentions=discord.AllowedMentions.none())
        await self.send_statistic(message, 'ban', member_obj, channel_obj, duration_display, updated_reason, expires_at, alias.alias_name, is_in_channel, bool(action), highest_role)
        
    # TODO
    
    async def handle_cow_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot cow the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}.')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to cow this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        select_sql = '''
            SELECT 1
            FROM active_cows
            WHERE discord_snowflake = $1
            AND channel_id = $2
        '''
        insert_cow_sql = '''
            INSERT INTO active_cows (guild_id, discord_snowflake, channel_id, created_at)
            VALUES ($1, $2, $3, $4)
            ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING
        '''
        insert_log_sql = '''
            INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
            VALUES ($1, $2, $3, $4, $5, $6)
            RETURNING created_at
        '''
        async with self.bot.db_pool.acquire() as conn:
            already_cowed = await conn.fetchval(select_sql, member_obj.id, channel_obj.id)
            if already_cowed:
                return await message.reply(content=f'\U0001F6AB {member_obj.mention} is already going vegan.', allowed_mentions=discord.AllowedMentions.none())
            created_at = await conn.fetchval(insert_log_sql, 'cow', member_obj.id, message.author.id, message.guild.id, channel_obj.id, 'Cowed a user')
            await conn.execute(insert_cow_sql, message.guild.id, member_obj.id, channel_obj.id, created_at)
            await message.reply(content=f'\U0001F525 {member_obj.mention} is going vegan!!! \U0001F525', allowed_mentions=discord.AllowedMentions.none())
    
    # TODO
    async def handle_flag_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None #'Tag a member or include their snowflake ID'
        prefix = args[1][0] if len(args) > 1 and args[1][:1] in {'+', '-', '='} else None
        reason = ' '.join(args[1:])[1:].lstrip() if prefix else ' '.join(args[1:]) if len(args) > 1 else 'No reason provided.'
 #'An optional reason. For modifications, `+` to append the reason, `=` to overwrite the reason, `-` to delete the reason')
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot flag the bot.')
        if member_obj:
            if member_obj.id in self.vegans:
                return await message.reply(content=f'\U0001F6AB You cannot flag a superhero.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        if not channel_obj:
            return await message.reply(content='\U0001F6AB Could not resolve a valid channel from the alias.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}.')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to flag this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        select_sql = '''
            SELECT reason
            FROM active_flags
            WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
        '''
        insert_sql = '''
            INSERT INTO active_flags (guild_id, discord_snowflake, channel_id, reason)
            VALUES ($1, $2, $3, $4)
            ON CONFLICT (guild_id, discord_snowflake, channel_id) DO NOTHING
        '''
        update_sql = '''
            UPDATE active_flags
            SET reason = $4
            WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
        '''
        is_coordinator = await is_owner_developer_coordinator_via_alias(message, alias)
        async with self.bot.db_pool.acquire() as conn:
            existing_flag = await conn.fetchrow(select_sql, message.guild.id, member_obj.id, channel_obj.id)
            reason_obj = Reason()
            if prefix in ['+', '=', '-']:
                if existing_flag:
                    reason_obj.load_old_reason(existing_flag['reason'])
                reason_obj.load_new_reason(reason)
                action = reason_obj.interpret_action(prefix)
                if action in ('delete', 'overwrite') and not is_coordinator:
                    return await message.reply(content='\U0001F6AB Only coordinators are allowed to overwrite the reason.')
                updated_reason = reason_obj.output_display()
            else:
                updated_reason = reason
            embed = discord.Embed(color=discord.Color.orange())
            embed.set_author(name=f'{member_obj.display_name} is flagged', icon_url=member_obj.display_avatar.url)
            embed.add_field(name='User', value=member_obj.mention, inline=True)
            embed.add_field(name='Channel', value=channel_obj.mention, inline=False)
            embed.add_field(name='Reason', value=update_reason, inline=False)
            await message.reply(embed=embed, allowed_mentions=discord.AllowedMentions.none())
   
    # TODO
    async def handle_role_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None #'Tag a member or include their snowflake ID'
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        is_owner_or_dev, is_coord_or_mod = await check_owner_dev_coord_mod(message, alias.channel.id)
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot give the bot a role.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        if not channel_obj:
            return await message.reply(content='\U0001F6AB Could not resolve a valid channel from the alias.')
        if not is_owner_or_dev and not is_coord_or_mod:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}')
        role_obj = message.guild.get_role(alias.role_id)
        if not role_obj:
            return await message.reply(content=f"\U000026A0\U0000FE0F Could not resolve role with ID `{alias.role_id}`.")
        if role_obj in member_obj.roles:
            return await message.reply(content=f'\U0001F6AB{member_obj.mention} already has {role_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        try:
            await member_obj.add_roles(role_obj, reason=f'Added role')
        except discord.Forbidden:
            return await message.reply(content=f'\U0001F6AB {member_obj.mention} was not successfully roled.', allowed_mentions=discord.AllowedMentions.none())
        return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} was given {role_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
    
    # TODO
    async def handle_text_mute_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None # 'Tag a member or include their snowflake ID'
        duration = args[1] if len(args) > 1 else '24h' # '(+|-)duration(m|h|d) \n 0 - permanent / 8h - default \n `+` to append, `-` to delete, `=` to overwrite reason'),
        reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.' #'Optional reason (required for 7 days or more)')
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj or member_obj.id in self.vegans:
            return await message.reply(content=f'\U0001F6AB Invalid target member: {member}.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB No permission to use `{alias.alias_name}` in {channel_obj.mention}.')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to mute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        is_coordinator = await is_owner_developer_coordinator_via_alias(message, alias)
        async with self.bot.db_pool.acquire() as conn:
            existing_text_mute = await conn.fetchrow('''
                SELECT expires_at, reason
                FROM active_text_mutes
                WHERE guild_id=$1 AND discord_snowflake=$2 AND channel_id=$3 AND room_name=$4
            ''', message.guild.id, member_obj.id, channel_obj.id, channel_obj.name)
            duration_obj = Duration()
            reason_obj = Reason()
            duration_obj.load_from_combined_duration_str(duration)
            if existing_text_mute and duration not in ['+', '=', '-']:
                expires_at = existing_text_mute['expires_at']
                old_reason = existing_text_mute['reason']
                reason_obj.load_old_reason(old_reason)
            if duration not in ['+', '=', '-']:
                expires_at = duration_obj.output_datetime()
                action = None
            else:
                reason_obj.load_new_reason(reason)
                action = reason_obj.interpret_action(duration_obj.get_prefix())
                if action in ('delete', 'overwrite') and not is_coordinator:
                    return await message.reply(content='\U0001F6AB Only coordinators are allowed to overwrite the reason.')
                updated_reason = reason_obj.output_display()
            if expires_at and expires_at <= now:
                return await message.reply(content='\U0001F6AB You cannot reduce a text-mute below the current time.')
            duration_obj.load_base(expires_at)
            duration_display = duration_obj.output_display()
            caps = await self.get_caps_for_channel(message.guild.id, channel_obj.id)
            active_cap = next((c for c in caps if c[0] == 'tmute'), None)
            if active_cap:
                cap_expires_at = active_cap[1]
            else:
                cap_expires_at = timedelta(days=7) + now
            if existing_text_mute and not expires_at < existing_text_mute['expires_at'] and not (is_coordinator or expires_at <= cap_expires_at):
                return await message.reply(content='\U0001F6AB Only coordinators can text-mute for longer than the channel cap.')
            is_in_channel = False
            if channel_obj:
                try:
                    await channel_obj.set_permissions(member_obj, send_messages=False, add_reactions=False)
                except discord.Forbidden:
                    return await message.reply(content=f'\U0001F6AB {member_obj.mention} was not successfully text-muted.', allowed_mentions=discord.AllowedMentions.none())
            await conn.execute('''
                INSERT INTO active_text_mutes (guild_id, discord_snowflake, channel_id, reason, expires_at, room_name)
                VALUES ($1,$2,$3,$4,$5,$6)
                ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name) DO UPDATE
                SET reason=$4, expires_at=$5
            ''', message.guild.id, member_obj.id, channel_obj.id, updated_reason, expires_at, channel_obj.name)
            await conn.execute('''
                INSERT INTO moderation_logs
                (action_type,target_discord_snowflake,executor_discord_snowflake,guild_id,channel_id,reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'textmute', member_obj.id, message.author.id, message.guild.id, channel_obj.id, f'Textmuted ({updated_reason})')
        embed = discord.Embed(
            title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} is text-muted",
            description=f"**User:** {member_obj.mention}\n**Channel/Room:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason}",
            color=discord.Color.orange()
        )
        await message.reply(embed=embed, allowed_mentions=discord.AllowedMentions.none())
        return await self.send_statistic(message, 'text_mute', member_obj, channel_obj, duration_display, updated_reason, expires_at, alias.alias_name, is_in_channel, bool(action), highest_role)
    
    # TODO
    async def handle_voice_mute_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None # 'Tag a member or include their snowflake ID'
        duration = args[1] if len(args) > 1 else '24h' # '(+|-)duration(m|h|d) \n 0 - permanent / 8h - default \n `+` to append, `-` to delete, `=` to overwrite reason'),
        reason = ' '.join(args[2:]) if len(args) > 2 else 'No reason provided.' #'Optional reason (required for 7 days or more)')
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        if not channel_obj:
            return await message.reply(content='\U0001F6AB Could not resolve a valid channel from the alias.')
        member_obj = await self.member_service.resolve_member(message, member)
        if member_obj:
            if member_obj.id in self.vegans:
                return await message.reply(content=f'\U0001F6AB You cannot mute a superhero.')
        if not member_obj:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot voice mute the bot.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB No permission to use `{alias.alias_name}` in {channel_obj.mention}.')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to mute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        is_coordinator = await is_owner_developer_coordinator_via_alias(message, alias)
        async with self.bot.db_pool.acquire() as conn:
            existing_mute = await conn.fetchrow('''
                SELECT expires_at, reason
                FROM active_voice_mutes
                WHERE guild_id = $1
                  AND discord_snowflake = $2
                  AND channel_id = $3
                  AND room_name = $4
            ''', message.guild.id, member_obj.id, static_channel_id, channel_obj.name)
            duration_obj = Duration()
            reason_obj = Reason()
            duration_obj.load_from_combined_duration_str(duration)
            if existing_mute and duration not in ['+', '=', '-']:
                expires_at = existing_mute['expires_at']
                old_reason = existing_mute['reason']
                reason_obj.load_old_reason(old_reason)
            if duration not in ['+', '=', '-']:
                expires_at = duration_obj.output_datetime()
                action = None
            else:
                reason_obj.load_new_reason(reason)
                action = reason_obj.interpret_action(duration_obj.get_prefix())
                if action in ('delete', 'overwrite') and not is_coordinator:
                    return await message.reply(content='\U0001F6AB Only coordinators are allowed to overwrite the reason.')
                updated_reason = reason_obj.output_display()
            if expires_at and expires_at <= now:
                return await message.reply(content='\U0001F6AB You cannot reduce a mute below the current time.')
            duration_obj.load_base(expires_at)
            duration_display = duration_obj.output_display()
            caps = await self.get_caps_for_channel(message.guild.id, channel_obj.id)
            active_cap = next((c for c in caps if c[0] == 'mute'), None)
            if active_cap:
                cap_expires_at = active_cap[1]
            else:
                cap_expires_at = timedelta(days=7) + now
            if existing_mute and not expires_at < existing_mute['expires_at'] and not (is_coordinator or expires_at <= cap_expires_at):
                return await message.reply(content='\U0001F6AB Only coordinators can mute for longer than the channel cap.')
        try:
            async with self.bot.db_pool.acquire() as conn:
                await conn.execute('''
                    INSERT INTO active_voice_mutes (guild_id, discord_snowflake, channel_id, expires_at, reason, target, room_name)
                    VALUES ($1, $2, $3, $4, $5, $6, $7)
                    ON CONFLICT (guild_id, discord_snowflake, channel_id, room_name, target)
                    DO UPDATE SET
                        expires_at = EXCLUDED.expires_at,
                        reason = EXCLUDED.reason
                ''', message.guild.id, member_obj.id, channel_obj.id, expires_at, updated_reason, 'user', channel_obj.name)
                await conn.execute('''
                    INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                    VALUES ($1, $2, $3, $4, $5, $6)
                ''', 'voice_mute', member_obj.id, message.author.id, message.guild.id, channel_obj.id, 'Voice muted a member')
        except Exception as e:
            logger.warning(f'DB insert failed: {e}')
            raise
        is_in_channel = False
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            is_in_channel = True
            try:
                await member_obj.edit(mute=True)
            except discord.Forbidden:
                return await message.reply(content=f"\U0001F6AB {member_obj.mention} was not successfully voice muted.", allowed_mentions=discord.AllowedMentions.none())
        embed = discord.Embed(
                title=f"{self.emoji.get_random_emoji()} {member_obj.display_name} is voice muted",
                description=f"**User:** {member_obj.mention}\n**Channel:** {channel_obj.mention}\n**Duration:** {duration_display}\n**Reason:** {updated_reason}",
                color=discord.Color.orange()
            )
        await message.reply(embed=embed, allowed_mentions=discord.AllowedMentions.none())
        await self.send_statistic(message, 'voice_mute', member_obj, channel_obj, duration_display, updated_reason, expires_at, alias.alias_name, is_in_channel, bool(action), highest_role)

    # TODO
    async def handle_unban_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot unban the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        if not channel_obj:
            return await message.reply(content='\U0001F6AB Could not resolve a valid channel from the alias.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to unban this `{highest_role}`.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT expires_at
                FROM active_bans
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND room_name = $4
            ''', message.guild.id, member_obj.id, channel_obj.id, channel_obj.name)
            if row and row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(message, alias):
                return await message.reply(content='\U0001F6AB Coordinator-only for undoing permanent bans.')
            try:
                if channel_obj:
                    await channel_obj.set_permissions(member_obj, overwrite=None)
            except discord.Forbidden:
                return await message.reply(content=f'\U0001F6AB {member_obj.mention} was not successfully unbanned.', allowed_mentions=discord.AllowedMentions.none())
            await conn.execute('''
                DELETE FROM active_bans
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND room_name = $4
            ''', message.guild.id, member_obj.id, channel_obj.id, channel_obj.name)
            await conn.execute('''
                INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                VALUES ($1,$2,$3,$4,$5,$6)
            ''', 'unban', member_obj.id, message.author.id, message.guild.id, channel_obj.id, f'Unbanned a user from room `{channel_obj.name}`')
        return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been unbanned from {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())

    # TODO
    async def handle_uncow_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj or not member:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot uncow the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to uncow this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        select_sql = '''
            SELECT 1
            FROM active_cows
            WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
        '''
        update_sql = '''
            DELETE FROM active_cows
            WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
        '''
        try:
            async with self.bot.db_pool.acquire() as conn:
                is_flagged = await conn.fetchval(select_sql, message.guild.id, member_obj.id, channel_obj.id)
                if not is_flagged:
                    return await message.reply(content=f'\U0001F6AB {member_obj.mention} has no active record in {channel_obj.mention}.')
                await conn.execute(update_sql, message.guild.id, member_obj.id, channel_obj.id)
                await conn.execute('''
                    INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                    VALUES ($1, $2, $3, $4, $5, $6)
                ''', 'uncow', member_obj.id, message.author.id, message.guild.id, channel_obj.id, 'Uncowed a user')
                await message.reply(content=f'👎 {member_obj.mention} is no longer going vegan. 👎', allowed_mentions=discord.AllowedMentions.none())
        except Exception as e:
            logger.warning(f'Database error occurred: {e}')
            raise
        
    # DONE
    async def handle_unflag_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj or not member:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot unflag the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias['channel_id'])
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to unflag this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        select_sql = '''
            SELECT 1
            FROM active_flags
            WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
        '''
        update_sql = '''
            DELETE FROM active_flags
            WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3
        '''
        try:
            async with self.bot.db_pool.acquire() as conn:
                is_flagged = await conn.fetchval(select_sql, message.guild.id, member_obj.id, channel_obj.id)
                if not is_flagged:
                    return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is not flagged for {channel_obj.mention}.')
                await conn.execute(update_sql, message.guild.id, member_obj.id, channel_obj.id)
                await conn.execute('''
                    INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason)
                    VALUES ($1, $2, $3, $4, $5, $6)
                ''', 'unflag', member_obj.id, message.author.id, message.guild.id, channel_obj.id, 'Unflagged a user')
                await message.reply(content=f'{self.emoji.get_random_emoji()} Unflagged {member_obj.mention} for channel {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        except Exception as e:
            logger.warning(f'Database error occurred: {e}')
            raise

    # TODO
    async def handle_unmute_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot unmute the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to unmute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT expires_at
                FROM active_voice_mutes
                WHERE guild_id = $1
                  AND discord_snowflake = $2
                  AND channel_id = $3
                  AND room_name = $4
                  AND target = 'user'
            ''', message.guild.id, member_obj.id, channel_obj.id, channel_obj.name)
            if not row:
                return await message.reply(content=f'\U0001F6AB {member_obj.mention} is not muted in {channel_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
            if row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(message, alias):
                return await message.reply(content='\U0001F6AB Coordinator-only for undoing permanent voice mutes.')
            await conn.execute('''
                DELETE FROM active_voice_mutes
                WHERE guild_id = $1
                  AND discord_snowflake = $2
                  AND channel_id = $3
                  AND room_name = $4
                  AND target = $5
            ''', message.guild.id, member_obj.id, channel_obj.id, channel_obj.name, 'user')
            if member_obj.voice and member_obj.voice.channel:
                try:
                    await member_obj.edit(mute=False)
                except discord.Forbidden:
                    return await message.reply(content=f'\U0001F6AB {member_obj.mention} was not successfully unmuted.', allowed_mentions=discord.AllowedMentions.none())
            await conn.execute('INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1,  $2, $3, $4, $5, $6)', 'unmute', member_obj.id, message.author.id, message.guild.id, channel_obj.id, 'Unmuted a member')
        if member_obj.voice and member_obj.voice.channel and member_obj.voice.channel.id == channel_obj.id:
            return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} has been unmuted in {channel_obj.mention}.',  allowed_mentions=discord.AllowedMentions.none())
        return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} is no longer marked as muted in {channel_obj.mention}.',  allowed_mentions=discord.AllowedMentions.none())

    # TODO
    async def handle_unrole_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
            return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot unrole the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias.channel.id)
        is_owner_or_dev, is_coord = await check_owner_dev_coord(message, channel_obj)
        if not is_coord and not is_owner_or_dev:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}.')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to derole this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        role_obj = message.guild.get_role(alias['role_id'])
        if not role_obj:
            return await message.reply(content=f"\U0001F6AB Could not resolve role with ID `{salias['role_id']}`.")
        if role_obj not in member_obj.roles:
            return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} does not have {role_obj.mention}.', allowed_mentions=discord.AllowedMentions.none())
        try:
            await member_obj.remove_roles(role_obj)
        except discord.Forbidden:
            return await message.reply(content=f'\U0001F6AB {member_obj.mention} was not successfully unroled.', allowed_mentions=discord.AllowedMentions.none())
        return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention} had {role_obj.mention} removed.', allowed_mentions=discord.AllowedMentions.none())
    
    # TODO
    async def handle_untextmute_alias(self, message: discord.Message, alias: Alias, args):
        await assert_owner_developer_coordinator_moderator(message, alias)
        member = args[0] if len(args) > 0 else None
        if not message.guild:
            return await message.reply(content='\U0001F6AB This command can only be used in servers.')
        member_obj = await self.member_service.resolve_member(message, member)
        if not member_obj:
             return await message.reply(content=f'\U0001F6AB Could not resolve a valid member from input: {member}.')
        if member_obj.bot:
            return await message.reply(content='\U0001F6AB You cannot undo a textmute on the bot.')
        channel_obj = await self.channel_service.resolve_channel(message, alias['channel_id'])
        if not channel_obj:
            return await message.reply(content='\U0001F6AB Could not resolve a valid channel from the alias.')
        is_owner_or_dev, is_mod_or_coord = await check_owner_dev_coord_mod(message, channel_obj)
        if not is_owner_or_dev and not is_mod_or_coord:
            return await message.reply(content=f'\U0001F6AB You do not have permission to use this command (`{alias.alias_name}`) in {channel_obj.mention}.')
        highest_role, success = await check_block(message, member_obj, channel_obj)
        if not success:
            return await message.reply(content=f'\U0001F6AB You are not allowed to un textmute this `{highest_role}` because they are a higher/or equivalent role than you in {channel_obj.mention}.')
        async with self.bot.db_pool.acquire() as conn:
            row = await conn.fetchrow('''
                SELECT expires_at
                FROM active_text_mutes
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND room_name = $4
            ''', message.guild.id, member_obj.id, channel_obj.id, channel_obj.name)
            if row and row['expires_at'] is None and not await is_owner_developer_coordinator_via_alias(message, alias):
                return await message.reply(content='\U0001F6AB Coordinator-only for undoing permanent text mutes.')
            try:
                await channel_obj.set_permissions(member_obj, send_messages=None)
            except discord.Forbidden:
                return await message.reply(content=f"\U0001F6AB {member_obj.mention} was not successfully untext-muted.", allowed_mentions=discord.AllowedMentions.none())
            await conn.execute('''
                DELETE FROM active_text_mutes
                WHERE guild_id = $1 AND discord_snowflake = $2 AND channel_id = $3 AND room_name = $4
            ''', message.guild.id, member_obj.id, channel_obj.id, channel_obj.name)
            await conn.execute('INSERT INTO moderation_logs (action_type, target_discord_snowflake, executor_discord_snowflake, guild_id, channel_id, reason) VALUES ($1, $2, $3, $4, $5, $6)', 'untmute', member_obj.id, message.author.id, message.guild.id, channel_obj.id, 'Untextmuted a user')
        return await message.reply(content=f'{self.emoji.get_random_emoji()} {member_obj.mention}\'s text muted in {channel_obj.mention} has been removed.', allowed_mentions=discord.AllowedMentions.none())

async def setup(bot: DiscordBot):
    cog = Aliases(bot)
    await bot.add_cog(cog)
