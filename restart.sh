set -e
cd ~/git/python/Vyrtuous
source ~/venv/bin/activate
poetry build --format wheel
pip uninstall vyrtuous
pip install dist/vyrtuous-6.0.6-py3-none-any.whl 
docker stop $(docker ps -aq) || true
docker rm $(docker ps -aq) || true
docker compose down -v --remove-orphans || true
docker container prune -f
docker volume prune -f
docker compose build
docker compose up -d
