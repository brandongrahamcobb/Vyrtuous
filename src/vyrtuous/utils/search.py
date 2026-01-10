import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.service.channel_service import ChannelService
from vyrtuous.service.check_service import role_check_without_specifics
from vyrtuous.service.member_service import MemberService
from vyrtuous.service.role_service import RoleService
from vyrtuous.enhanced_member.enhanced_member import EnhancedMember
from vyrtuous.moderation_action.moderation_action import ModerationAction
from vyrtuous.moderation_action.voice_mute import VoiceMute
from vyrtuous.room.room import Room

def generate_skipped_set_pages(chunk_size, field_count, pages, skipped, title):
    embed = discord.Embed(
        title=title,
        description='\u200b',
        color=discord.Color.blue()
    )
    lines = []
    for snowflake in skipped:
        if field_count >= chunk_size:
            embed.description = '\n'.join(lines)
            pages.append(embed)
            embed = discord.Embed(
                title=f'{title} continued...',
                color=discord.Color.red()
            )
            lines = []
            field_count = 0
        lines.append(str(snowflake))
        field_count += 1
    embed.description = '\n'.join(lines)
    pages.append(embed)
    return pages

def generate_skipped_dict_pages(chunk_size, field_count, pages, skipped, title):
    for guild_snowflake, list in skipped.items():
        embed = discord.Embed(
            color=discord.Color.red(),
            title=f'{title} ({guild_snowflake})'
        )
        field_count = 0
        lines = []
        for snowflake in list:
            if field_count >= chunk_size:
                embed.description = '\n'.join(lines)
                pages.append(embed)
                embed = discord.Embed(
                    color=discord.Color.red(),
                    title=f'{title} ({guild_snowflake}) continued...'
                )
                field_count = 0
                lines = []
            lines.append(str(snowflake))
            field_count += 1
        embed.description = '\n'.join(lines)
        pages.append(embed)

async def send_pages(obj, pages, state):
    if pages:
        try:
            return await state.end(success=pages)
        except Exception as e:
            try:
                return await state.end(warning=f'\U000026A0\U0000FE0F Embed size is too large. Limit the scope.')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
    else:
        try:
            return await state.end(warning=f'\U000026A0\U0000FE0F No {obj.PLURAL.lower()} found.')
        except Exception as e:
            return await state.end(error=f'\u274C {str(e).capitalize()}')

class Search:

    def __init__(self):
        self.channel_service = ChannelService()
        self.member_service = MemberService()
        self.role_service = RoleService()

    async def match_scope_to_objects(self, channel_obj, guild_obj, member_obj, role_obj, scope, obj):
        match scope:
            case 'all':
                objs = await obj.fetch_all()
            case 'channel':
                if not channel_obj or not guild_obj:
                    raise 
                objs = await obj.fetch_by_channel_and_guild(
                    channel_snowflake=channel_obj.id,
                    guild_snowflake=guild_obj.id
                )
            case 'guild':
                if not guild_obj:
                    raise
                objs = await obj.fetch_by_guild(
                    guild_snowflake=guild_obj.id
                )
            case 'member':
                if not guild_obj or not member_obj:
                    raise
                objs = await obj.fetch_by_guild_and_member(
                    guild_snowflake=guild_obj.id,
                    member_snowflake=member_obj.id
                )
            case 'role':
                if not guild_obj or not role_obj:
                    raise
                objs = await obj.fetch_by_guild_and_role(
                    guild_snowflake=guild_obj.id,
                    role_snowflake=role_obj.id
                )
            case _:
                if not channel_obj or not guild_obj:
                    objs = await obj.fetch_by_channel_and_guild(
                        channel_snowflake=channel_obj.id,
                        guild_snowflake=guild_obj.id
                    )
        return objs
    
    @classmethod
    async def resolve_objs_from_target(cls, ctx_interaction_or_message, obj, state, target):
        bot = DiscordBot.get_instance()
        channel_obj, guild_obj, member_obj, role_obj = None, None, None, None
        scope = None
        highest_role = await role_check_without_specifics(
            ctx_interaction_or_message=ctx_interaction_or_message
        )
        if target and isinstance(target, str) and target.lower() == 'all':
            if highest_role not in ('System Owner', 'Developer'):
                try:
                    return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                        f'You are not authorized to list {obj.PLURAL.lower()} ' \
                        f'across all servers.'
                    )
                except Exception as e:
                    return await state.end(error=f'\u274C {str(e).capitalize()}')
            scope = 'all'
        elif target:
            try:
                channel_obj = await cls.channel_service.resolve_channel(
                    ctx_interaction_or_message=ctx_interaction_or_message,
                    channel_str=target
                )
                scope = 'channel'
            except Exception as e:
                try:
                    member_obj = await cls.member_service.resolve_member(
                        ctx_interaction_or_message=ctx_interaction_or_message,
                        member_str=target
                    )
                    scope = 'member'
                except Exception as e:
                    guild_obj = bot.get_guild(int(target))
                    if not guild_obj:
                        try:
                            role_obj = await cls.role_service.resolve_role(
                                ctx_interaction_or_message=ctx_interaction_or_message,
                                role_str=target
                            )
                            scope = 'role'
                        except Exception as e:
                            try:
                                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                                    f'Scope must be one of: `all`, channel ID/mention, ' \
                                    f'member ID/mention, server ID or empty. Received: {target}.'
                                )
                            except Exception as e:
                                return await state.end(error=f'\u274C {str(e).capitalize()}')
                    else:
                        if highest_role in (
                            'Coordinator',
                            'Moderator',
                            'Everyone'
                        ):
                            try:
                                return await state.end(warning=f'\U000026A0\U0000FE0F ' \
                                    f'You are not authorized to list {obj.PLURAL.lower()} ' \
                                    f'for specific servers.'
                                )
                            except Exception as e:
                                return await state.end(error=f'\u274C {str(e).capitalize()}')
                        scope = 'guild'
        else:
            scope = None
            channel_obj = ctx_interaction_or_message.channel
            guild_obj = ctx_interaction_or_message.guild
            if isinstance(ctx_interaction_or_message.guild, (commands.Context, discord.Message)):
                member_obj = ctx_interaction_or_message.author
            else:
                member_obj = ctx_interaction_or_message.user
    
        match obj:
            case EnhancedMember():
                if obj() == Administrator():
                    match scope:
                        case 'all':
                            objs = await obj.fetch_all()
                        case 'guild':
                            if not guild_obj:
                                raise
                            objs = await obj.fetch_by_guild(
                                guild_snowflake=guild_obj.id
                            )
                        case 'member':
                            if not guild_obj or not member_obj:
                                raise
                            objs = await obj.fetch_by_guild_and_member(
                                guild_snowflake=guild_obj.id,
                                member_snowflake=member_obj.id
                            )
                        case 'role':
                            if not guild_obj or not role_obj:
                                raise
                            objs = await obj.fetch_by_guild_and_role(
                                guild_snowflake=guild_obj.id,
                                role_snowflake=role_obj.id
                            )
                        case _:
                            if not guild_obj:
                                objs = await obj.fetch_by_guild(
                                    guild_snowflake=guild_obj.id
                                )
                            msg = f'No {obj.PLURAL.lower()} exist for ' \
                                f'{guild_obj.name}'
                else:
                    objs = await cls.match_scope_to_objects(
                        channel_obj=channel_obj,
                        guild_obj=guild_obj,
                        member_obj=member_obj,
                        role_obj=role_obj,
                        obj=obj,
                        scope=scope
                    )
                    msg = f'No {obj.PLURAL.lower()} exist for ' \
                        f'{channel_obj.mention} in {guild_obj.name}'
            case ModerationAction():
                if obj() == VoiceMute():
                    target = 'user'
                    match scope:
                        case 'all':
                            objs = await obj.fetch_all_by_target(target=target)
                        case 'channel':
                            objs = await obj.fetch_by_channel_guild_and_target(
                                channel_snowflake=channel_obj.id,
                                guild_snowflake=ctx_interaction_or_message.guild.id,
                                target=target
                            )
                        case 'guild':
                            objs = await obj.fetch_by_guild_and_target(
                                guild_snowflake=guild_obj.id,
                                target=target
                            )
                        case 'member':
                            objs = await obj.fetch_by_guild_member_and_target(
                                guild_snowflake=ctx_interaction_or_message.guild.id,
                                member_snowflake=member_obj.id,
                                target=target
                            )
                        case _:
                            objs = await obj.fetch_by_channel_guild_and_target(
                                channel_snowflake=ctx_interaction_or_message.channel.id,
                                guild_snowflake=ctx_interaction_or_message.guild.id,
                                target=target
                            )
                            msg = f'No {obj.PLURAL.lower()} exist for ' \
                                f'{channel_obj.mention} in {guild_obj.name}'
                else:
                    objs = await cls.match_scope_to_objects(
                        channel_obj=channel_obj,
                        guild_obj=guild_obj,
                        member_obj=member_obj,
                        role_obj=role_obj,
                        obj=obj,
                        scope=scope
                    )
                    msg = f'No {obj.PLURAL.lower()} exist for ' \
                        f'{channel_obj.mention} in {guild_obj.name}'
            case Room():
                objs = await cls.match_scope_to_objects(
                        channel_obj=channel_obj,
                        guild_obj=guild_obj,
                        member_obj=member_obj,
                        role_obj=role_obj,
                        obj=obj,
                        scope=scope
                    )
                msg = f'No {obj.PLURAL.lower()} exist for ' \
                        f'{channel_obj.mention} in {guild_obj.name}'
            case _:
                return
        
        if not isinstance(objs, list):
            objs = [objs]
        if not objs:
            try:
                if scope:
                    msg = f'No {obj.PLURAL.lower()} exist for scope: {scope}.'
                return await state.end(warning=f'\U000026A0\U0000FE0F {msg}')
            except Exception as e:
                return await state.end(error=f'\u274C {str(e).capitalize()}')
        return objs
        
            

