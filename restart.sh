#!/bin/bash
set -e
cd ~/git/python/Vyrtuous
source ~/venv/bin/activate
poetry build --format wheel
pip uninstall vyrtuous
pip install dist/vyrtuous-6.0.6-py3-none-any.whl 
docker stop $(docker ps -aq) || true
docker compose build
docker compose up -d db
docker exec -i $(docker ps -qf "name=db") psql -U postgres -d vyrtuous < ./schema/vyrtuous.sql
docker compose run --service-ports vyrtuous
