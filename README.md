**The Vyrtuous Project** is a vegan-owned Discord bot written in Python. It brings together:

* PostgreSQL database information about mutes.
* Loading and unloading Cogs, the fundamental components of the Discord bot.

## Features

Moderation Features

• `Room Muting` — Mutes in Discord rooms only persist in the room itself, not the whole server.


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

## Commands

Replace `<prefix>` with your configured prefix (default `?`).

• `?help <command>` or `?help`

• `?load <path-to-cog>`

• `?reload <path-to-cog>`

• `?sync <~|^|*|>` or `?sync`

## License

Distributed under GPL-3.0-or-later. See LICENSE for details.
