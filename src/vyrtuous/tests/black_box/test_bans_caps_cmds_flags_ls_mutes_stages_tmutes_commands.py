"""test_chown_temp_temps_xtemp_commands.py The purpose of this program is to black box test the temporary room commands.
Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from typing import Optional
from vyrtuous.inc.helpers import ROLE_ID
from vyrtuous.tests.black_box.test_suite import (
    bot,
    config,
    extract_embed_text,
    guild,
    not_privileged_author,
    prepared_command_handling,
    prefix,
    privileged_author,
    text_channel,
    voice_channel_one,
    RESET, YELLOW, RED, GREEN
)
from vyrtuous.utils.emojis import EMOJIS
import pytest


def parse_duration_seconds(token, unit_map):
    magnitude = token[1:-1]
    unit = token[-1]
    if not magnitude.isdigit():
        return None
    seconds = int(magnitude)
    factor = {"s": 1, "m": 60, "h": 3600, "d": 86400}.get(unit_map.get(unit))
    return seconds * factor if factor else None


def extract_duration_seconds(cmd, unit_map):
    for part in cmd.split():
        if part and part[0] in "=+-":
            seconds = parse_duration_seconds(part, unit_map)
            if seconds is not None:
                return seconds
    return None


def apply_active_cap(alias_cases, cap_seconds, unit_map):
    out = []
    for cmd, permission, ref_channel, ref_guild, ref_member, should_warn in alias_cases:
        warn = should_warn
        duration = extract_duration_seconds(cmd, unit_map)
        if duration is not None and duration > cap_seconds:
            warn = True
        out.append((cmd, permission, ref_channel, ref_guild, ref_member, warn))
    return out


def generate_cap_commands(durations):
    moderation_types = ["ban", "vmute", "tmute"]
    caps = []
    for mod_type in moderation_types:
        for duration in durations:
            caps.append(
                (
                    "Administrator",
                    f"cap {{voice_channel_one_id}} {mod_type} {duration}",
                    True,
                    False,
                    False,
                    False,
                )
            )
    return caps


NORMAL_ALIAS_SETUP = "alias {alias_type} {alias} {voice_channel_one_id}"
ROLE_ALIAS_SETUP = "alias {alias_type} {alias} {voice_channel_one_id} {role_id}"
ALIAS_CREATION_VARIANTS = NORMAL_ALIAS_SETUP
ALIAS_TEARDOWN = ["xalias {alias}"]
DAY_UNITS = {"d", "day", "days"}
DURATIONS = {"0", "1", "24", "48"}
HOUR_UNITS = {"h", "hr", "hrs", "hour", "hours"}
LIST_ALIASES = {
    # "bans": [
    #     dict(
    #         undo=False,
    #         alias="testban",
    #         alias_type="ban",
    #         all_list_permission_role="Developer",
    #         list_permission_role="Administrator",
    #         alias_permission_role="Moderator",
    #         timed=True,
    #     ),
    #     dict(
    #         undo=True,
    #         alias="testunban",
    #         alias_type="unban",
    #         all_list_permission_role="Developer",
    #         list_permission_role="Administrator",
    #         alias_permission_role="Moderator",
    #         timed=False,
    #     ),
    # ],
    # "flags": [
    #     dict(
    #         undo=False,
    #         alias="testflag",
    #         alias_type="flag",
    #         all_list_permission_role="Developer",
    #         list_permission_role="Administrator",
    #         alias_permission_role="Moderator",
    #         timed=False,
    #     ),
    #     dict(
    #         undo=True,
    #         alias="testunflag",
    #         alias_type="unflag",
    #         all_list_permission_role="Developer",
    #         list_permission_role="Administrator",
    #         alias_permission_role="Moderator",
    #         timed=False,
    #     ),
    # ],
    "ls": [
        dict(
            undo=False,
            alias="testv",
            alias_type="vegan",
            all_list_permission_role="Developer",
            list_permission_role="Administrator",
            alias_permission_role="Moderator",
            timed=False,
        ),
        dict(
            undo=True,
            alias="testc",
            alias_type="carnist",
            all_list_permission_role="Developer",
            list_permission_role="Administrator",
            alias_permission_role="Moderator",
            timed=False,
        ),
    ],
    # "mutes": [
    #     dict(
    #         undo=False,
    #         alias="testvm",
    #         alias_type="vmute",
    #         all_list_permission_role="Developer",
    #         list_permission_role="Administrator",
    #         alias_permission_role="Moderator",
    #         timed=True,
    #     ),
    #     dict(
    #         undo=True,
    #         alias="testunvm",
    #         alias_type="unvmute",
    #         all_list_permission_role="Developer",
    #         list_permission_role="Administrator",
    #         alias_permission_role="Moderator",
    #     ),
    # ],
    # "tmutes": [
    #     dict(
    #         undo=False,
    #         alias="testtm",
    #         alias_type="tmute",
    #         all_list_permission_role="Developer",
    #         list_permission_role="Administrator",
    #         alias_permission_role="Moderator",
    #         timed=True,
    #     ),
    #     dict(
    #         undo=True,
    #         alias="testuntm",
    #         alias_type="untmute",
    #         all_list_permission_role="Developer",
    #         list_permission_role="Administrator",
    #         alias_permission_role="Moderator",
    #         timed=False,
    #     ),
    # ],
}
LIST_CASES = {
    # "cmds": dict(
    #     command=None,
    #     list_command="cmds",
    #     all_list_permission_role="Developer",
    #     list_permission_role="Moderator",
    #     command_permission_role=None,
    # ),
    # "stages": dict(
    #     command="stage",
    #     list_command="stages",
    #     all_list_permission_role="Developer",
    #     list_permission_role="Moderator",
    #     command_permission_role="Administrator",
    # ),
    # "temps": dict(
    #     command="temp",
    #     list_command="temps",
    #     all_list_permission_role="Developer",
    #     list_permission_role="Moderator",
    #     command_permission_role="Administrator",
    # ),
    # "vrs": dict(
    #     command="vr",
    #     list_command="vrs",
    #     alias_type="vrs",
    #     all_list_permission_role="Developer",
    #     list_permission_role="Moderator",
    #     command_permission_role="Administrator",
    # ),
}
MINUTE_UNITS = {"m", "min", "mins", "minute", "minutes"}
PREFIXES = {"", "=", "+", "-"}
REASON = "test"
ROLE_ALIASES = [
    dict(
        alias="testrole",
        alias_type="role",
        all_list_permission_role="Developer",
        list_permission_role="Administrator",
        alias_permission_role="Coordinator",
    ),
    dict(
        alias="testunrole",
        alias_type="unrole",
        all_list_permission_role="Developer",
        list_permission_role="Administrator",
        alias_permission_role="Coordinator",
    ),
]
ROLE_ALIAS_SETUP = ROLE_ALIAS_SETUP
SCOPE = ["all", "{voice_channel_one_mention}", "{voice_channel_one_id}", "{guild_id}"]
SECOND_UNITS = {"s", "sec", "secs", "second", "seconds"}
UNIT_MAP = {
    **dict.fromkeys(DAY_UNITS, "d"),
    **dict.fromkeys(HOUR_UNITS, "h"),
    **dict.fromkeys(MINUTE_UNITS, "m"),
    **dict.fromkeys(SECOND_UNITS, "s"),
}

CAP_SECONDS = parse_duration_seconds("=24h", UNIT_MAP)


def build_alias_cases(
    LIST_ALIASES,
    LIST_CASES,
    PREFIXES, 
    DURATIONS,
    UNIT_MAP,
    REASON,
    SCOPE,
    NORMAL_ALIAS_SETUP,
):
    cases = []
    durations = sorted(DURATIONS, key=lambda x: int(x))
    scopes = sorted(SCOPE)
    units = sorted(UNIT_MAP.keys())

    def add_case(
        cmd, permission, ref_channel, ref_guild, ref_member, should_warn, phase
    ):
        cases.append(
            {
                "cmd": cmd,
                "permission": permission,
                "ref_channel": ref_channel,
                "ref_guild": ref_guild,
                "ref_member": ref_member,
                "should_warn": should_warn,
                "phase": phase,
            }
        )

    def setup_alias(alias):
        cmd = NORMAL_ALIAS_SETUP.format(
            alias_type=alias["alias_type"],
            alias=alias["alias"],
            voice_channel_one_id="{voice_channel_one_id}",
        )
        add_case(
            cmd, alias["alias_permission_role"], True, False, False, False, "setup"
        )

    def apply_alias(alias, undo_alias):
        if not alias.get("timed"):
            cmd = f"{alias['alias']} {{member_id}} {REASON}"
            add_case(
                cmd, alias["alias_permission_role"], True, False, True, False, "apply"
            )
            return
        prefixes = ["=", "+", "-"]
        for u in units:
            for d in durations:
                for p in prefixes:
                    if d == "0" and p == "-":
                        continue
                    cmd = f"{alias['alias']} {{member_id}} {p}{d}{u} {REASON}"
                    add_case(
                        cmd,
                        alias["alias_permission_role"],
                        True,
                        False,
                        True,
                        False,
                        "apply",
                    )
                    if d == "0":
                        undo_cmd = f"{undo_alias['alias']} {{member_id}}"
                        add_case(
                            undo_cmd,
                            undo_alias["alias_permission_role"],
                            True,
                            False,
                            True,
                            False,
                            "undo",
                        )

    def list_permutations(list_name, alias_type, permission, should_warn):
        for scope in scopes:
            ref_guild = scope in ("all", "{guild_id}")
            add_case(
                f"{list_name} {scope}",
                permission,
                False,
                ref_guild,
                False,
                should_warn,
                "list",
            )

    def undo_alias(alias):
        cmd = f"{alias['alias']} {{member_id}}"
        add_case(cmd, alias["alias_permission_role"], True, False, True, False, "undo")

    for list_name, aliases in LIST_ALIASES.items():
        normal = next(a for a in aliases if not a["undo"])
        undo = next(a for a in aliases if a["undo"])
        setup_alias(normal)
        setup_alias(undo)
        apply_alias(normal, undo)
        permission = (
            LIST_CASES[list_name]["list_permission_role"]
            if list_name in LIST_CASES
            else normal["list_permission_role"]
        )
        list_permutations(list_name, normal["alias_type"], permission, True)

        undo_alias(undo)
        permission = (
            LIST_CASES[list_name]["list_permission_role"]
            if list_name in LIST_CASES
            else normal["list_permission_role"]
        )
        list_permutations(list_name, normal["alias_type"], permission, True)

    return cases


# ALIAS_CASES=build_alias_cases(LIST_ALIASES, LIST_CASES, PREFIXES, DURATIONS, DAY_UNITS, HOUR_UNITS, MINUTE_UNITS, SECOND_UNITS, UNIT_MAP, REASON, SCOPE, NORMAL_ALIAS_SETUP)
ALIAS_CASES = build_alias_cases(
    LIST_ALIASES,
    LIST_CASES,
    PREFIXES,
    DURATIONS,
    UNIT_MAP,
    REASON,
    SCOPE,
    NORMAL_ALIAS_SETUP,
)

ALIAS_CASES_TUPLE = []

for case in ALIAS_CASES:
    ALIAS_CASES_TUPLE.append(
        (
            case["cmd"],
            case["permission"],
            case["ref_channel"],
            case["ref_guild"],
            case["ref_member"],
            case["should_warn"],
        )
    )

FINAL_CASES = []
FINAL_CASES.extend(ALIAS_CASES_TUPLE)
FINAL_CASES.extend(generate_cap_commands(DURATIONS))
FINAL_CASES.extend(apply_active_cap(ALIAS_CASES_TUPLE, CAP_SECONDS, UNIT_MAP))


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "cmd,permission,ref_channel,ref_guild,ref_member,should_warn", FINAL_CASES
)
async def test_bans_caps_cmds_flags_ls_mutes_stages_tmutes_commands(
    bot,
    cmd: Optional[str],
    guild,
    not_privileged_author,
    permission,
    prefix: Optional[str],
    privileged_author,
    ref_channel,
    ref_guild,
    ref_member,
    should_warn,
    text_channel,
    voice_channel_one,
):
    channel_values = (voice_channel_one.mention, voice_channel_one.id)
    guild_values = (guild.name, guild.id)
    member_values = (not_privileged_author.mention, not_privileged_author.id)
    formatted = cmd.format(
        voice_channel_one_id=voice_channel_one.id,
        voice_channel_one_mention=voice_channel_one.mention,
        member_id=not_privileged_author.id,
        guild_id=guild.id,
        role_id=ROLE_ID,
    )
    captured = await prepared_command_handling(
        author=privileged_author,
        bot=bot,
        channel=text_channel,
        content=formatted,
        guild=guild,
        highest_role=permission,
        prefix=prefix,
    )
    message = captured[0]["message"]
    if message.embeds:
        embed = message.embeds[0]
        content = extract_embed_text(embed)
    elif message.embed:
        content = extract_embed_text(message.embed)
    else:
        content = message.content
    message_type = captured[0]["type"]
    if message_type == "error":
        print(f"{RED}Error:{RESET} {content}")
    if message_type == "warning":
        print(f"{YELLOW}Warning:{RESET} {content}")
        if should_warn:
            assert True
        else:
            assert False
    if message_type == "success":
        print(f"{GREEN}Success:{RESET} {content}")
        if ref_channel:
            assert any(
                str(channel_value) in content for channel_value in channel_values
            )
        if ref_guild:
            assert any(str(guild_value) in content for guild_value in guild_values)
        if ref_member:
            assert any(str(member_value) in content for member_value in member_values)
        # if moderation_type:
        #     assert moderation_type in content
        assert any(emoji in content for emoji in EMOJIS)
