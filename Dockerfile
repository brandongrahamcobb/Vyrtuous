FROM python:3.12-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    gnupg \
    lsb-release && \
    echo "deb http://apt.postgresql.org/pub/repos/apt $(lsb_release -cs)-pgdg main" > /etc/apt/sources.list.d/pgdg.list && \
    wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | apt-key add - && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    python3-dev \
    libglib2.0-dev \
    libgl1-mesa-glx \
    qtbase5-dev \
    libqt5gui5 \
    libqt5webkit5-dev \
    libqt5svg5-dev \
    postgresql-client-17 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    POETRY_VIRTUALENVS_CREATE=false

# Set working directory
WORKDIR /app
COPY pyproject.toml poetry.lock* ./
# Install Poetry
RUN pip install --no-cache-dir poetry

RUN poetry install --no-root


CMD ["python", "-u", "vyrtuous/main.py"] 
