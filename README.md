**The Vyrtuous Project** is a vegan-owned Discord bot written in Python. It brings together:

* Loading and unloading Cogs, the fundamental components of the Discord bot.
* Permission scaling for server developers, channel coordinators and channel moderators.
* PostgreSQL database information about ban, flags and mutes.

## Features

Moderation Features

• `Room Ban` — Bans in Discord rooms only persist in the room itself, not the whole server.
It includes ban durations including a permanent ban.
It also includes an optional reason.

• `Room Muting` — Mutes in Discord rooms only persist in the room itself, not the whole server.
It also includes an optional reason.

* `alias <alias_type> <alias_name> <channel>` - Creates a ban, flag, mute, unban or unmute alias for a specific channel
* `aliases <channel>` - Lists all the aliases present in a channel.
* `coord <channel> <member>`
* `coords` - Lists coordinators in the server where the command was run.
* `dev <member>` - Creates a developer in the server where the command was run.
* `devs` - Lists developers in the server where the command was run.
* `flags <channel>` - Lists members who are flagged in the channel.
* `help <command>` - Interactive help command paginating commands for members.
* `mod <channel> <member>` - Creates a moderator in the channel specified.
* `mods <channel>` - Lists all mods in the channel specified.
* `roleban <channel> <role>` - Syncs a channel's ban alias with an assigned moderated role.'
* `xalias <alias_name>` - Deletes an alias.
* `xcoord <channel> <member>` - Deletes a coordinator from a specified channel.
* `xdev <member>` - Deletes a developer in the server it was run.
* `xmod <channel> <member>` - Deletes a moderator in a channel.

Lifecycle Features

• `load <path-to-cog>` — Loads the bot's cogs after an unload.

• `reload <path-to-cog>` - Reloads the bot's cogs on the fly.

• `sync <~|^|*|>` or `sync` - Syncs the command tree for the bot to a guild for application command access.

## Installation

Prerequisites:

• Python 3.13+ (`python3 --version`)

• pip (`pip --version`)

• PostgreSQL (`psql --version`)

• Docker (`docker --version`)

1. Clone the repo

```bash
git clone https://github.com/brandongrahamcobb/Vyrtuous.git
```

3. Install python

4. Install pip

5. Create & activate a virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
```

6. Install Docker

7. Run restart.sh
```bash
./.restart.sh
```
8. Subsequent restarts should run start.sh
Run restart.sh
```bash
./.start.sh
```

## Configuration

On first run, The Vyrtuous Project will prompt you to enter and confirm:

• Discord bot token

```txt
How many api keys do you want to make?
1
"Discord"
Enter api_key
```

• Discord command prefix
• Discord character limit (regular = 2000, nitro = 4000)
• Discord owner ID
• Discord test guild

Your settings are saved this in the container:

```txt
~/.config/vyrtuous/config.yml
```

Subsequent launches read from this file—no environment variables needed.

The bot will load or create its config, connect to Discord, and register commands.

## License

Distributed under GPL-3.0-or-later. See LICENSE for details.
