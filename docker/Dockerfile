FROM --platform=linux/x86_64 mambaorg/micromamba:1.4.5 as app

ARG MENTALIST_VERSION='1.0.0'

USER root
WORKDIR /

# Install requirements
RUN apt-get update && \
    apt-get install -y \
    git \
    build-essential \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Create environment
RUN micromamba install --yes --name base \
    --channel conda-forge \
    --channel bioconda \
    biopython \
    julia=1.9.1 \
    hdf5 \
    mpfr \
    libxml2 \
    unzip && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN julia -e 'using Pkg; Pkg.add([ "Distributed", "ArgParse", "JSON", "DataStructures", "JLD", "GZip", "Blosc", "FileIO", "TextWrap", "LightXML", "JuMP", "Gurobi"])' && \
    julia -e 'using Pkg; Pkg.add(name="BioSequences", version="1.0")'

RUN git clone https://github.com/WGS-TB/MentaLiST.git && \
    cd /MentaLiST && \
    git checkout docker -- && \
    chmod +x /MentaLiST/src/mentalist 

WORKDIR /data

ENV PATH="${PATH}:/opt/conda/bin:/MentaLiST/src:/MentaLiST/scripts"
