# MentaLiST

MentaLiST -- The MLST pipeline developed by the PathOGiST research group

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

## Quick Start

A notebook with basic commands to install MLST schema and running the allele calling algorithm can be found [here](docs/Basic%20Usage.ipynb).
