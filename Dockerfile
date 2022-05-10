FROM python:3.9

# set environment varibles
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

WORKDIR /app
ENV PATH="/app:${PATH}"

COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . ./

