**The Vyrtuous Project** is a vegan-owned Discord bot written in Python. It brings together:

* Loading and unloading Cogs, the fundamental components of the Discord bot.
* Permission scaling for server mute admins, server developers, channel coordinators and channel moderators.
* PostgreSQL database information about ban, flags and mutes.
* Aliasing to enable room by room moderation.
* Reengineered server mute as local mute.
* Duration and expiration timers with reasons.

## Project Flow
![Vyrtuous](resources/images/VyrtuousUML.svg)

## Features

Moderation Features

• `Room Ban` — Bans in Discord rooms only persist in the room itself, not the whole server.
It includes ban durations including a permanent ban.
Coordinators and above can use the permanent ban.
It also includes an optional reason (non-optional if permanent).

• `Room Muting` — Mutes in Discord rooms only persist in the room itself, not the whole server.
It includes mute durations including a permanent mute.
Coordinators and above can use the permanent ban.
It also includes an optional reason.

• `Text Muting` — Mutes in Discord text-channels only persist in the room itself, not the whole server.
It includes mute durations including a permanent mute.
Coordinators and above can use the permanent ban.
It also includes an optional reason.

## Static Commands
* `admin <member>` - Creates a server muter. Requires an owner.
* `alias <alias_type> <alias_name> <channel>` - Creates a ban, cow, flag, mute, unban , uncow or unmute alias for a specific channel.
* `cmds <channel>` - Lists all the aliases present in a channel or use `all` to list all channels in a guild.
* `coord <channel> <member>` - Creates a coordinator for a given channel. Permitted to use by developers.
* `coords` - Lists coordinators in the server where the command was run or use `all` to list all coordinators in a guild.
* `dev <member>` - Creates a developer in the server where the command was run. Requires bot owner or server owner permission.
* `devs` - Lists developers in the server where the command was run.
* `flags <channel>` - Lists members who are flagged in the channel or use `all` to list all flags in a guild.
* `help <command>` - Interactive help command paginating commands for members.
* `mod <channel> <member>` - Creates a moderator in the channel specified requires a developer or above.
* `mods <channel>` - Lists all mods in the channel specified or use `all` to list all moderators in a guild.
* `xalias <alias_name>` - Deletes an alias. Requires a coordinator.
* `xadmin <member>` - Deletes a server muter. Requires an owner.
* `xcoord <channel> <member>` - Deletes a coordinator from a specified channel. Requires a developer or above.
* `xdev <member>` - Deletes a developer in the server it was run.
* `xmod <channel> <member>` - Deletes a moderator in a channel.

## Alias Commands
* `ban <member> <duration> <reason>` - Bans a member for a certain timeframe (24h by default) via 30m, 2h, 1d like strings for duration and an optional reason for moderaotrs. Permanent bans are set using 0 for duration and it requires a coordinator.
* `flag <member> <reason>` - Flags a member. Requires a moderator.
* `mute <member> <duration> <reason>` - Voice mutes a member for a certain timeframe (24h by default) via 30m, 2h, 1d like strings for duration and an optional reason for moderaotrs. Permanent voice mutes are set using 0 for duration and it requires a coordinator.
* `tmute <member> <duration> <reason>` - Text mutes a member for a certain timeframe (24h by default) via 30m, 2h, 1d like strings for duration and an optional reason for moderaotrs. Permanent text-mutes are set using 0 for duration and it requires a coordinator.
* `unban <member>    - Unbans a member for a room. Requires moderator.
* `unflag <member>   - Unflags a member for a room. Requires moderator. 
* `unmute <member>   - Unvoicemutes a member for a room. Requires moderator. 
* `untmute <member>  - Untextmutes a member for a room. Requires moderator. 

Lifecycle Features

• `backup` — Creates an immediate backup and sends it over Discord.

• `load <path-to-cog>` — Loads the bot's cogs after an unload.

• `reload <path-to-cog>` - Reloads the bot's cogs on the fly.

• `sync <~|^|*|>` or `sync` - Syncs the command tree for the bot to a guild for application command access.

## Installation

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
4. Profit

The bot will connect to Discord, and register commands.

## License

Distributed under GPL-3.0-or-later. See LICENSE for details.
