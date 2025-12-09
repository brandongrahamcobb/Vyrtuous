FROM python:3.12-alpine

RUN apk add --no-cache postgresql-client

ENV POETRY_VIRTUALENVS_CREATE=false

WORKDIR /app
RUN mkdir -p /app/backups

RUN pip install --upgrade pip
RUN pip install --no-cache-dir poetry

COPY pyproject.toml poetry.lock* README.md ./
RUN poetry install

COPY healthcheck.sh ./
COPY src ./

CMD ["python", "-u", "vyrtuous/main.py"]
