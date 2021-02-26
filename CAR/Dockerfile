FROM frolvlad/alpine-miniconda3

WORKDIR /code

# Install conda environment from yaml
COPY environment.yml environment.yml
RUN conda env create -f environment.yml && \
    conda clean -afy && \ 
    find /opt/conda/ -follow -type f -name '*.a' -delete &&\ 
    find /opt/conda/ -follow -type f -name '*.pyc' -delete &&\
    find /opt/conda/ -follow -type f -name '*.js.map' -delete
        
# Get packages for rdkit
RUN apk update && \
    apk --no-cache add bash ca-certificates wget libxrender libxext nodejs npm && \
    update-ca-certificates && \
    npm i npm@latest -g 

# Install the node stuff from package.json
# NB also issue with react-beatiful-dnd and dependency with old
# ver of react therfor used '--legacy-peer-deps' tag
COPY package.json .
RUN npm install --quiet --legacy-peer-deps 

# Copy over all the files
COPY . . 

