**The Vyrtuous Project** is a vegan-owned Discord bot written in Python. It brings together:

![Vyrtuous UML Diagram](resources/pictures/VyrtuousUML.svg)

* Custom room mute commands.
* Custom room unmute commands.
* No conflict with manual muting.
* PostgreSQL database information like developers, moderators and roles.
* Loading and unloading Cogs, the fundamental components of the Discord bot.
* Backup redelployment commands.

## Features

Moderation Features

• `alias <mute|unmute> <alias> <channelId>` — Create a room specific alias for users to use to mute or unmute members in their room.

• `dev <userId> <guildId>` — Give a user permissions to operate the bot in the guild as a developer.

• `devs` — Lists all the current devs in the guild where the command was run.

• `help <command>` - Get command-specific usage information.

• `mod <userId> <channelId>` — Give a user's permissions in a room to mute and unmute members.

• `mods` — Lists all the current room mods in the guild where the command was run.

• `xalias <mute|unmute> <alias>` — Delete a room specific alias for mute or unmute.

• `xdev <userId> <guildid>` — Revoke a user's permissions to operate the bot in the guild as a developer.

• `xmod <userId> <channelId>` — Revoke a user's permissions in a room to mute and unmute members.


Lifecycle Features

• `load <path-to-cog>` — Loads the bot's cogs after an unload.

• `reload <path-to-cog>` - Reloads the bot's cogs on the fly.

• `sync <~|^|*|>` or `sync` - Syncs the command tree for the bot to a guild for application command access.

## Installation

Prerequisites:
• Python 3.13+ (`python3 --version`)
• pip (`pip --version`)
• PostgreSQL (`psql --version`)

1. Clone the repo

```bash
git clone https://github.com/brandongrahamcobb/Vyrtuous.git
```

2. Create & activate a virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
```

3. Inside the venv, install Poetry

```bash
pip install poetry
```

4. Build the wheel

```bash
poetry build --format wheel
pip install dist/vyrtuous-0.0.1-py3-none-any.whl
```

5. Create the PostgreSQL database

```bash
createdb vyrtuous
```

6. Run the SQL setup script

```bash
psql -U postgres -d vyrtuous -f vyrtuous.sql
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

Your settings are saved to:

```txt
~/.config/vyrtuous/config.yml
```

Subsequent launches read from this file—no environment variables needed.

## Running

With your venv active, simply run:

```bash
vyrtuous
```

The bot will load or create its config, connect to Discord, and register commands.

## Commands

Replace `<prefix>` with your configured prefix (default `?`).

• `?alias <mute|unmute> <alias> <channelId>`

• `?dev <userId>`

• `?devs`

• `?help <command>` or `?help`

• `?load <path-to-cog>`

• `?mod <userId> <channelId>`

• `?mods`

• `?reload <path-to-cog>`

• `?xalias <mute|unmute> <alias>`

• `?xdev <userId>`

• `?xmod <userId> <channelId>`

• `?sync <~|^|*|>` or `?sync`

## License

Distributed under GPL-3.0-or-later. See LICENSE for details.
