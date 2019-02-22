#!/bin/sh

cp -r $SRC_DIR/src/*.jl $PREFIX/bin
cp -r $SRC_DIR/scripts $PREFIX
ln -s $PREFIX/bin/MentaLiST.jl $PREFIX/bin/mentalist
chmod +x $PREFIX/bin/mentalist

julia -e 'using Pkg; Pkg.add([ "Distributed", "ArgParse", "BioSequences", "JSON", "DataStructures", "JLD", "GZip", "Blosc", "FileIO", "TextWrap", "LightXML", "JuMP", "Gurobi"])'

rm -f "$PREFIX"/share/julia/site/lib/v*/*.ji
rm -rf "$PREFIX"/share/julia/site/v*/METADATA
rm -f "$PREFIX"/share/julia/site/v*/META_BRANCH
