#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(RNAmrf))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))


filemrf=args[1]
fileseq=args[2]
fileout=args[3]

wt_h=as.numeric(args[4])
wt_j=as.numeric(args[5])

gap_e=as.numeric(args[6])
gap_o=as.numeric(args[7])

# seqfiles=list.files(fileseqs,pattern = "raw")

seqs=seqinr::read.fasta(fileseq)
mymrf=RNAmrf::read_mrf_renum(filemrf)


  
a2b_1 = align_seq2mrf(seqs[[1]], mymrf, wt_h=wt_h,wt_j=wt_j,gap_ext=gap_e,gap_open=gap_o, debug = FALSE)
a2b_2 = align_seq2mrf(seqs[[2]], mymrf, wt_j=wt_j,wt_h=wt_h,gap_ext=gap_e,gap_open=gap_o, debug = FALSE)
seqaln = pair_a2b2aln(a2b_1, a2b_2, seqs = seqs)

seqinr::write.fasta(seqaln,
                  names = names(seqs),
                  file.out = file.path(fileout, paste0(basename(fileseq), ".mrfaln.fa")))

