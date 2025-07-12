**The Vyrtuous Project** is a vegan-owned Discord bot written in Python. It brings together:

![Vyrtuous UML Diagram](resources/pictures/VyrtuousUML.svg)

* Room moderator data persistence
* Custom room mute commands
* Custom room unmute commands
* Loading and unloading Cogs, the fundamental components of the Discord bot
* Backing up roles and redeploying the backup

## Features

Moderation Features

• `delalias <mute|unmute> <alias>` — Delete a room specific alias for mute or unmute.

• `setalias <mute|unmute> <alias> <channelId>` — Create a room specific alias for users to use to mute or unmute members in their room.

• `help <command>` - Get command-specific usage information.

• `give_mod <userId> <channelId>` — Give a user's permissions in a room to mute and unmute members.

• `list_mods <userId>` — Lists all the current room mods in the server where the command was run.

• `revoke_mod <userId> <channelId>` — Revoke a user's permissions in a room to mute and unmute members.

• `give_dev <userId>` — Give a user permissions to operate the bot in the server as a developer.

• `list_devs <userId>` — Lists all the current devs in the server where the command was run.

• `revoke_dev <userId>` — Revoke a user's permissions to operate the bot in the server as a developer.

Lifecycle Features
• `load <path-to-cog>` — Loads the bot's cogs after an unload.

• `reload <path-to-cog>` - Reloads the bot's cogs on the fly.

• `sync <~|^|*|>` or `sync` - Syncs the command tree for the bot to a server for application command access.

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
/venv/lib/<python-version>/site-packages/vyrtuous/.config/config.yml
```

Subsequent launches read from this file—no environment variables needed.

## Running

With your venv active, simply run:

```bash
vyrtuous
```

The bot will load or create its config, connect to Discord, and register commands.

## Commands

Replace `<prefix>` with your configured prefix (default `!`).

• `!delalias <mute|unmute> <alias>`
• `!setalias <mute|unmute> <alias> <channelId>`
• `!help <command>`
• `!give_mod <userId> <channelId>`
• `!list_mods <userId>`
• `!revoke_mod <userId> <channelId>`
• `!give_dev <userId>`
• `!list_devs <userId>`
• `!revoke_dev <userId>`
• `!load <path-to-cog>`
• `!reload <path-to-cog>`
• `!sync <~|^|*|>` or `sync`

## License

Distributed under GPL-3.0-or-later. See LICENSE for details.
