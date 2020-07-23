FROM julia:0.6.3

RUN apt-get update && apt-get install -y libhdf5-dev

ADD https://github.com/WGS-TB/MentaLiST/archive/v0.2.4.tar.gz /tmp/mentalist.tar.gz

RUN mkdir -p /opt/mentalist /tmp/mentalist &&\
tar -xvf /tmp/mentalist.tar.gz --strip-components=1 -C /tmp/mentalist &&\
mv /tmp/mentalist/src /opt/mentalist/ &&\
#mv /tmp/mentalist/*.toml /opt/mentalist/ &&\
mv /tmp/mentalist/scripts/* /opt/mentalist/ &&\
rm -r /tmp/mentalist* &&\
echo '#!/usr/bin/env bash \n\
exec julia "/opt/mentalist/src/MentaLiST.jl" "$@" \n' > "/usr/bin/mentalist" &&\
chmod +x "/usr/bin/mentalist"

ENV JULIA_DEPOT_PATH="/usr/share/julia/site"

RUN julia -e '\
Pkg.update();\
Pkg.add("Bio");\
Pkg.add("OpenGene");\
Pkg.add("Logging");\
Pkg.add("ArgParse");\
Pkg.add("Lumberjack");\
Pkg.add("FastaIO");\
Pkg.add("JLD");\
Pkg.add("DataStructures");'