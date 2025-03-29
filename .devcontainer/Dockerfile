FROM debian:bullseye-slim

# Copy Python requirements and NPM packages
COPY ./src/requirements-dev.txt /container/requirements-dev.txt

USER root

WORKDIR /container

# Update and install
RUN \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get install --no-install-recommends -y gnupg2 python3-rdkit python3-pip python3-pandas python3-psycopg2 nodejs npm git git-crypt nano locales build-essential && \ 
 apt-get autoclean -y && \
 pip install -r requirements-dev.txt  

# # Set the locale to the UK
RUN locale-gen en_GB.UTF-8 && \
    update-locale LANG en_GB.UTF-8

ENV LANG=en_GB.UTF-8
ENV LANGUAGE=en_GB:en
