FROM python:3.12-slim

# RUN apt-get update && \
#     apt-get install -y --no-install-recommends \
#         wget \
#         gnupg \
#         lsb-release \
#         ca-certificates && \
#     mkdir -p /etc/apt/keyrings && \
#     wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc \
#         | gpg --dearmor -o /etc/apt/keyrings/postgres.gpg && \
#     echo "deb [signed-by=/etc/apt/keyrings/postgres.gpg] http://apt.postgresql.org/pub/repos/apt $(lsb_release -cs)-pgdg main" \
#         > /etc/apt/sources.list.d/pgdg.list && \
#     apt-get update && \
#     apt-get install -y --no-install-recommends postgresql-client-18 && \
#     apt-get clean && rm -rf /var/lib/apt/lists/*

ENV POETRY_VIRTUALENVS_CREATE=false

WORKDIR /app
RUN mkdir -p /app/backups

COPY src README.md pyproject.toml poetry.lock* ./

RUN pip install --upgrade pip
RUN pip install --no-cache-dir poetry

RUN which python

RUN poetry install
CMD ["python", "-u", "vyrtuous/main.py"]
