FROM debian:bullseye-slim

# Copy over all the files
COPY requirements.txt /app/

WORKDIR /app

# Update and install
RUN \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get install --no-install-recommends -y python3-rdkit python3-pip python3-pandas python3-psycopg2 nodejs npm && \ 
 apt-get autoclean -y &&\
 pip install -r requirements.txt &&\  
 mkdir logs &&\
 cd logs &&\
 touch logfile.log

WORKDIR /app/CAR/frontend

RUN \
 npm i npm@latest -g && \
 npm install --quiet --legacy-peer-deps
