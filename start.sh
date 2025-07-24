#!/bin/bash
docker stop $(docker ps -aq) || true
docker compose up -d
docker exec -i $(docker ps -qf "name=db") psql -U postgres -d vyrtuous < ./schema/password_protected/encrypted.sql
docker compose run --service-ports vyrtuous
