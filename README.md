# MentaLiST

MentaLiST -- The MLST pipeline developed by the PathOGiST research group

## Build status

[![Build Status](https://travis-ci.org/WGS-TB/MentaLiST.svg?branch=mentalist_v0.1)](https://travis-ci.org/WGS-TB/MentaLiST)

## Abstract

MLST (multi-locus sequence typing) is a classic technique for genotyping
bacteria, widely applied for pathogen outbreak surveillance. Traditionally,
MLST is based on identifying strain types from a set of a small number of
housekeeping genes. With the increased availability of whole-genome sequencing
(WGS) data, MLST methods have evolved toward larger typing schemes, based on a
few hundred genes (core genome MLST, cgMLST) to a few thousand genes (whole
genome MLST, wgMLST). Such large-scale MLST schemes have been shown to provide
a finer resolution and are increasingly used in various contexts such as
hospital outbreaks or foodborne pathogen outbreaks. This methodological shift
raises new computational challenges, especially given the large size of the
schemes involved. Very few available MLST callers are currently capable of
dealing with large cgMLST and wgMLST schemes.

We introduce MentaLiST, a new MLST caller, based on a k-mer counting algorithm and written in the Julia language, specifically designed and implemented to handle large typing schemes. Tests on real and simulated data to show that MentaLiST is faster than any other available MLST caller while providing the same or better accuracy, and is capable of dealing with MLST schema with up to thousands of genes while requiring limited computational resources. 

## Installation

### Linux

The easiest way of installing MentaLiST is by creating a new environment with [Conda](https://conda.io/docs/). To create a new conda environment that includes MentaLiST, run:
```
conda create -n mentalist -c bioconda mentalist
```

then activate it by running:
```
source activate mentalist
```

Once the mentalist conda environment is active, you should be able to run MentaLiST. Typing:
```
mentalist -h 
```
produces the help output. 

The conda environment can be deactivated by running:
```
 source deactivate 
```
### macOS (+10.8)

There is currently no conda recipe for julia on macOS, so a more manual installation process is required.

1. Download the [julia-0.5.2.dmg](https://julialang-s3.julialang.org/bin/mac/x64/0.5/julia-0.5.2-mac64.dmg) file from julialang.org and install it by dragging the `Julia-0.5.app` bundle into your Applications folder.
2. Launch the `Julia-0.5` application and install all of the dependencies listed in the [REQUIRE](REQUIRE) file:

```julia
julia> Pkg.update()
julia> Pkg.add("Bio")
julia> Pkg.add("OpenGene")
julia> Pkg.add("Logging")
julia> Pkg.add("ArgParse")
julia> Pkg.add("Lumberjack")
julia> Pkg.add("Suppressor")
```
3. Add `/Applications/Julia-0.5.app/Contents/Resources/julia/bin/julia` to your `PATH` 
4. Clone the MentaLiST git repostory (https://github.com/WGS-TB/MentaLiST.git). MentaLiST can be run directly from the repository: `src/mentalist -h`.

## Quick Start

A notebook with basic commands to install MLST schema and running the allele calling algorithm can be found [here](docs/Basic%20Usage.ipynb).

You can also read more detailed instruction on how to deal with [novel alleles](docs/Novel%20allele%20detection%20with%20MentaLiST.ipynb).
