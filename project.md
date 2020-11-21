
====================================================
## MRFalign project

### Done
- Bralibase 2.1 K2, 27 families without pseudoknot, pairwise alignment test
    - Results 
        - MRFalign can generate good alignment for low sequence ID, but slightly worse than CMalign

## Issues
- How should I prove that MRFalign is better than CMalign at aligning RNA with pk? 
    - RFAM Data may be biased: only conserved PKs are recorded.
        - RFAM only record RNAs with conserved "PK", so CM can already align those PKs well based on sequence similarity
        - The trained MRF model may not have "enough" data
     
### Possible plans

- Find better test set for pseudoknot RNA alignment
- 

- Synthetic MSA experiment 
    1. Given RNA SS: generate MSA(train/test) that satisfies SS constraints: basepairs in canonical, loop regions can be randomly mutated.  
    2. Construct MRF and CM from MSA, align test sequence with MRF and CM model
    3. Benchmark alignment with reference

- ncRNA homolog search with CMalign and MRF
    - Reduce the search space to smaller genomes ()

====================================================
Thesis writing plan

- Chapter1 Intro

- Chapter2 DFIRE_RNA 

- Chapter3 RNAcmap

- Chapter4 RNA SS prediction with DCA constraints

- Chapter5 RNA MRF alignment

- Chapter6 Summary