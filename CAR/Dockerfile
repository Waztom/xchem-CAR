FROM continuumio/miniconda3

WORKDIR /code

# Try get npm on system?
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y nodejs npm

# Install conda environment from yaml
ADD environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml

# Need to update npm otherwise the build fails at 'npm install' 
# due to compatibility issues with version of node js
RUN npm i npm@latest -g   

# Install the node stuff from package.json
# NB also issue with react-beatiful-dnd and dependency with old
# ver of react therfor used '--legacy-peer-deps' tag
COPY package.json .
RUN npm install --quiet --legacy-peer-deps 

# Copy over all the files
COPY . . 