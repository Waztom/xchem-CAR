FROM debian:bullseye-slim

# Copy Python and Node reqs
COPY ./CAR/package.json /workspace/package.json
COPY ./CAR/requirements-dev.txt /workspace/requirements-dev.txt

# Update args in docker-compose.yaml to set the UID/GID of the "vscode" user.
ARG USER_UID=1000
ARG USER_GID=$USER_UID
RUN if [ "$USER_GID" != "1000" ] || [ "$USER_UID" != "1000" ]; then groupmod --gid $USER_GID vscode && usermod --uid $USER_UID --gid $USER_GID vscode; fi

WORKDIR /workspace

# Update and install
RUN \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get install --no-install-recommends -y gnupg2 python3-rdkit python3-pip python3-pandas python3-psycopg2 nodejs npm git git-crypt nano && \ 
 apt-get autoclean -y &&\
 pip install -r requirements-dev.txt &&\  
 npm i npm@latest -g && \
 npm install --quiet --legacy-peer-deps 

