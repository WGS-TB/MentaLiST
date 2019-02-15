# MentaLiST

MentaLiST -- The MLST pipeline developed by the PathOGiST research group

## Build status

[![Build Status](https://travis-ci.org/WGS-TB/MentaLiST.svg?branch=mentalist_v0.2)](https://travis-ci.org/WGS-TB/MentaLiST)
[![Coverage Status](https://coveralls.io/repos/github/WGS-TB/MentaLiST/badge.svg?branch=mentalist_v0.2)](https://coveralls.io/github/WGS-TB/MentaLiST?branch=mentalist_v0.2)

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

There are two ways to install MentaLiST. The simplest way is with `conda`, which installs Julia and Python environments with all the necessary requirements automatically. Currently, this is only available for Linux. For macOS, please follow the manual installation. 

The other way is by manually downloading and installing Julia, and the required libraries.

### Conda installation

**Note: We've received reports that the current conda builds of MentaLiST are broken. MentaLiST was recently updated to run on Julia 1.1, so these builds are quite out-of-date in respect with the current GitHub version. Until this is resolved, we recommend installing MentaLiST with the manual installation, as described below. **

The easiest way of installing MentaLiST is by creating a new environment with [Conda](https://conda.io/docs/). 

Ensure that you have both the `conda-forge` and `bioconda` channels enabled:

```
conda config --show channels
channels:
  - conda-forge
  - bioconda
  - defaults
```

To create a new conda environment called `mentalist` that includes MentaLiST, run:

```
conda create -n mentalist mentalist
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

### Manual Installation

MentaLiST is written in Julia, and currently runs on versions >= 1.0. If you don't have it already, install an [official Julia-1.1 binary](https://julialang.org/downloads/) or more recent version.

Then, clone this repo to obtain the MentaLiST source. After that, run Julia-1.1 and install the required packages listed in the REQUIRE file. For macOS, you'll have to drag the `Julia-1.1.app` bundle from the .dmg file into your Applications folder, and Launch the app.
You can also add `/Applications/Julia-1.1.app/Contents/Resources/julia/bin/julia` to your `PATH` 

```
julia
julia> using Pkg
julia> Pkg.add("Distributed")
julia> Pkg.add("ArgParse")
julia> Pkg.add("BioSequences")
julia> Pkg.add("JSON")
julia> Pkg.add("DataStructures")
julia> Pkg.add("JLD")
julia> Pkg.add("GZip")
julia> Pkg.add("Blosc")
julia> Pkg.add("FileIO")
julia> Pkg.add("TextWrap")
julia> Pkg.add("LightXML")
```
#### Optional Packages
If you want to use the allele coverage method (`-c` option), useful for compressing very large wgMLST schemes when building a new MLST database, you will need a copy of Gurobi ILP solver, and also the Gurobi and JuMP julia packages. Please check https://github.com/JuliaOpt/Gurobi.jl for installation instructions. 

After installing the dependencies, `MentaLiST` can be run directly from the repository: `src/mentalist -h`.

## Quick Start

A notebook with basic commands to install MLST schema and running the allele calling algorithm can be found [here](docs/Basic%20Usage.ipynb).

You can also read more detailed instruction on how to deal with [novel alleles](docs/Novel%20allele%20detection%20with%20MentaLiST.ipynb).
