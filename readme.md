
## Dataset Description 

This repo contains the data for paper RNAmrf and scripts to reproduce figures in the paper.

Structure of the dir (Use RF00165 for example)

- `data/RFAM_PK`
    - `data/RFAM_PK/RFAM/RF00165`: an RFAM family
        - input:
            - `input/RF00165.a2m`: reference MSA in [a2m format](http://emboss.sourceforge.net/docs/themes/AlignFormats.html): aligned positions are in uppercase and "-", insertions are lowercase and ".". converted from RF00165.sto by [esl-reformat](http://eddylab.org/infernal/)
            - `input/RF00165.afa`: reference MSA in [afa format](https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/multalignviewer/afasta.html): a fasta-like format.
            - `input/RF00165.afa.mrf`: MRF model produced by GREMLIN
            - `input/RF00165.unaligned.fasta`: unaligned sequences
        - output:
            - `cmalign/RF00165.cmalign.a2m`: cmalign aligned MSA in a2m format
            - `rnamrf/RF00165.rnamrf.a2m`: rnamrf aligned MSA in a2m format
        
- `data/synthetic` follows similar directory structures.

Other files:

- apsi.txt : average pairwise sequence index, can be calculated by `R/calc_apsi.R`
- CLEN.txt : length of the reference alignment
- neff.txt : number of effective sequence, calculated using [GREMLIN](https://github.com/sokrypton/GREMLIN_CPP)
- sci.scores: structural conservation index, calculated using [RNAalifold](https://www.tbi.univie.ac.at/RNA/)


## How to Use RNAmrf

[RNAmrf](https://github.com/tcgriffith/RNAmrf/) is an open-source (MIT) R package to align RNA sequences to a Markov Random Field model. You can install the development version from github:

```
## install from GitHub
remotes::install_github('tcgriffith/RNAmrf')
```

The script 'R/RUN_mrfaln.R' provides a wrapper to the RNAmrf, which can be executed from commandline:

```
## In directory ./example
## Usage: RUN_mrfaln.R input.mrf input.fasta output.a2m
../R/RUN_mrfaln.R input/example.afa.mrf input/unaligned.fasata rnamrf/example.rnamrf.a2m
```

The alignment quality can be benchmarked:

```bash
## In the main directory
./R/bench_pkaln_dir.R ./example

## The benchmark compares the alignment rnamrf/example.rnamrf.a2m with reference MSA input/example.a2m. 

## output:
         col_all  col_loop  pair_all pair_pk pair_nonpk method  caseid
rnamrf 0.8300427 0.7533383 0.8541667  0.8475  0.8608333 rnamrf example

#col_all:      alignment accuracy of all columns (percentage of correctly aligned positions comparing test and reference MSA)
#col_loop:     alignment accuracy of non-paired columns
#pair_all:     alignment accuracy of base-pairs
#pair_pk :     alignment accuracy of pseudoknots
#pair_nonpk :  alignment accuracy of non-pseudoknot stem
```

This shows that RNAmrf has an 0.83 alignment accuracy, both Pseudoknot and non-pseudoknot stems have an average alignment accuracy of 0.85.
