FROM julia:1.9.1

RUN apt-get update && \
    apt-get install -y git && \
    rm -rf /var/lib/apt/lists/*

RUN julia -e 'using Pkg; Pkg.add([ "Distributed", "ArgParse", "BioSequences", "JSON", "DataStructures", "JLD", "GZip", "Blosc", "FileIO", "TextWrap", "LightXML", "JuMP", "Gurobi"])'

RUN git clone https://github.com/WGS-TB/MentaLiST.git && \
    cd /MentaLiST && \
    git checkout docker && \
    chmod +x /MentaLiST/src/mentalist 

ENV PATH="${PATH}:/MentaLiST/src"