**The Vyrtuous Project** is a vegan-owned Discord bot written in Python. It brings together:

* Permission scaling for guild owners, guild administrators, channel coordinators and channel moderators.
* Aliasing to enable channel by channel moderation.
* Channel-scoped right-click mutes.
* Temporary channels.
* Visibility of all moderation actions.
* Video-only channels.

# Project Flow
![Vyrtuous](resources/images/VyrtuousUML.svg)

# Features

### Aliases
* `ban <member> <duration> <reason>` - Bans a member.
* `flag <member> <reason>` - Flags a member.
* `vegan <member>` - Tracks new vegans.
* `tmute <member> <duration> <reason>` - Text mutes a member.
* `vmute <member> <duration> <reason>` - Voice mutes a member.
* `carnist <member>` - Removes a vegan.
* `unban <member>` - Unbans a member.
* `unflag <member>` - Unflags a member.
* `untmute <member>` - Undoes a text-mutes for a member.
* `unvmute <member>` - Undoes a voice-mute for a member.

## Static Commands

## Guild Owner
* `arole <role>` - Promotes all members part of a role to administrator along with active tracking of who retains the role or loses the role, revoking or granting administrator if a role removed or given, respectively.
* `hero <member>` - Grants or revokes invincibility for a member.

## Administrator
* `alias <alias_type> <alias_name> <channel>` - Creates a ban, cow, flag, mute, unban , uncow or unmute alias for a specific channel.
* `aroles` - Lists all current administrator roles.
* `cap <channel> <'ban', 'tmute' or 'vmute`> <hours> - Limits the duration of moderator actions by moderators (Coordinators and above bypass these caps).
* `chown <channel> <member>` - Transfers temporary room ownership to a new member.
* `clear <member>` - Removes all moderation actions active on a member.
* `coord <channel> <member>` - Creates a coordinator for a given channel.
* `mtrack <channel> <scope> <entry_type> <snowflakes>` - Setup for moderation action tracking by 'general', 'channel' or 'member' and their respective ID's.
* `rmv <source_channel> <target_channel>` - Moves all members from one voice channel to another voice channel.
* `smute <member>` - Toggles a server mute for a member.
* `temp <channel> <member>` - Toggles a temporary channel with a channel owner for maintaining temporary channel aliases and their respective bans, flags, text-mutes, vegans and voice-mutes.
* `temps <scope>` - Lists temporary channels and their aliases.
* `vr <channel>` - Toggles a video-only channel for maintaining video-only chaennl aliases and their respective bans, flags, text-mutes, vegans and voice-mutes.
* `vrs <scope>` - Lists video-only channels and their aliases.
* `stage <channel> <duration>` - Toggles a psuedo-stage for a voice channel.
* `track <channel>` - Toggles moderation action tracking for a given channel.
* `tracks <channel>` - Lists the setup created by `!mtrack`.
* `xalias <alias_name>` - Deletes an alias.

## Coordinator
* `mod <channel> <member>` - Toggles a moderator in the channel specified.
* `rmute <channel>` - Mutes everyone in a room besides yourself.
* `xrmute <channel>` - Unmutes everyone in a room.

## Moderator
* `caps <scope>` - Lists moderation duration caps.
* `cmds <channel>` - Lists aliases.
* `del <message> <channel> ` - Deletes a message.
* `flags <scope>` - Lists active flags.
* `migrate <channel_name> <channel> ` - Transfers all the previous setup for a temporary channel to a new channel.
* `mstage <member> <channel>` - Toggles a voice-mute for a member in a stage channel.
* `mutes <scope>` - Lists active voice-mutes.
* `stages <scope>` - Lists active stages.
* `summary <member> <scope>` - Reports all moderation actions on a member.
* `tmutes <scope>` - Lists active text-mutes.

## Everyone
* `admins <scope>` - Lists administrators.
* `coords <scope>` - Lists coordinators.
* `help <command>` - Interactive help command paginating commands for members.
* `ls <scope>` - Lists vegans.
* `mods <scope>` - Lists moderators.
* `roleid <role>` - Gets the role ID.
* `survey <channel>` - Gets all the elevated members in a channel.

# Installation

Prerequisites:

- Python 3.13+ (`python3 --version`)
- Poetry (https://python-poetry.org/docs/#installation)
- Docker (`docker --version`)

1. Clone and cd into directory

```bash
git clone https://gitlab.com/vyrtuous/vyrtuous
cd vyrtuous
```

2. Copy `.env.example` to `.env` and populate the variables
```bash
cp .env.example .env
nano .env
```

3. Use compose to build and run the stack
```sh
docker compose up
```

The bot will connect to Discord, and register commands.

# License

Distributed under GPL-3.0-or-later. See LICENSE for details.
