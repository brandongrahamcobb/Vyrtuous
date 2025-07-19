#!/bin/bash

set -e
cd ~/git/python/Vyrtuous

# 1. Activate virtual environment
source ~/venv/bin/activate

docker stop $(docker ps -aq) || true
docker compose build
docker compose up -d db
docker exec -i $(docker ps -qf "name=db") psql -U postgres -d vyrtuous < ./schema/vyrtuous.sql
docker compose run --service-ports vyrtuous
