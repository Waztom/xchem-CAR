FROM debian:bullseye-slim

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
 mkdir logs &&\
 cd logs &&\
 touch logfile.log

# Expose port
EXPOSE 8000

# Run on Entry
CMD cd /app/ ; npm run build ; python3 manage.py wait-for-db ; python3 manage.py migrate backend  ; python3 manage.py makemigrations backend ; python3 manage.py runserver 0.0.0.0:8000