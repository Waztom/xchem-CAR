FROM debian:bullseye-slim

# Copy Python requirements and NPM packages
COPY ./src/requirements-dev.txt /container/requirements-dev.txt

USER root

WORKDIR /container

# Install base dependencies and NVM
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install --no-install-recommends -y \
    curl \
    gnupg2 \
    procps \
    python3-rdkit \
    python3-pip \
    python3-psycopg2 \
    git \
    git-crypt \
    nano \
    locales \
    build-essential && \
    apt-get autoclean -y && \
    pip install -r requirements-dev.txt

# Install NVM and Node.js 18
ENV NVM_DIR=/root/.nvm
ENV NODE_VERSION=18.19.0

RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.7/install.sh | bash && \
    . $NVM_DIR/nvm.sh && \
    nvm install $NODE_VERSION && \
    nvm alias default $NODE_VERSION && \
    nvm use default

# Set the locale to the UK
RUN locale-gen en_GB.UTF-8 && \
    update-locale LANG en_GB.UTF-8

ENV LANG=en_GB.UTF-8
ENV LANGUAGE=en_GB:en
