FROM debian:bullseye-slim

# Copy Python and Node reqs
COPY package.json /app/package.json
COPY requirements.txt /app/requirements.txt

WORKDIR /app

# Copy over all the files
COPY . . 

# Update and install
RUN \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get install --no-install-recommends -y python3-rdkit python3-pip python3-pandas python3-psycopg2 nodejs npm && \ 
 apt-get autoclean -y &&\
 pip install -r requirements.txt &&\  
 npm i npm@latest -g && \
 npm install --quiet --legacy-peer-deps &&\
 cd /app &&\
 npm run build

# Expose port
EXPOSE 8000

# Run on Entry
CMD cd /app/CAR ; python ./manage.py migrate ; python ./manage.py runserver 0.0.0.0:8000