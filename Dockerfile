FROM python:3.9-slim

RUN apt-get -qq update \
    && apt-get install -y --no-install-recommends \
       procps \
    && rm -rf /var/lib/apt/lists/*

# set environment varibles
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

WORKDIR /app
ENV PATH="/app:${PATH}"

COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . ./

