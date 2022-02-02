FROM debian:bullseye-slim

# Copy Python and Node reqs
COPY ./CAR/frontend/package.json /workspace/CAR/frontend/package.json
COPY ./CAR/requirements-dev.txt /workspace/requirements-dev.txt

USER root

WORKDIR /workspace

# Update and install
RUN \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get install --no-install-recommends -y gnupg2 python3-rdkit python3-pip python3-pandas python3-psycopg2 nodejs npm git git-crypt nano && \ 
 apt-get autoclean -y &&\
 pip install -r requirements-dev.txt  

# Install frontend npm packages
WORKDIR /workspace/CAR/frontend

RUN \
 npm i npm@latest -g && \
 npm install --quiet --legacy-peer-deps

CMD chown -R root /workspace/CAR/frontend