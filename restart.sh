#!/bin/bash

set -e
cd /Users/spawd/git/python/Vyrtuous

# 1. Activate virtual environment
source /Users/spawd/venv/bin/activate

# 2. Bump version
VERSION=$(python src/vyrtuous/utils/inc/bump_version.py | tail -n 1)

# 3. Export version so shell has access
export VYRTUOUS_VERSION="$VERSION"

# 4. Build and install wheel
python -m poetry build --format wheel
pip install --force-reinstall "dist/vyrtuous-$VYRTUOUS_VERSION-py3-none-any.whl"

# 5. Drop/recreate schema
psql -U postgres -f main.sql
psql -U postgres -d vyrtuous -f vyrtuous.sql

# 6. Load updated vars into current shell
source /Users/spawd/.bashrc_vyrtuous

