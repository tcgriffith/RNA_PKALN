#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(RNAmrf))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))


filemrf=args[1]
fileseq=args[2]
# fileout=args[3]
fileout="./"

# wt_h=as.numeric(args[4])
# wt_j=as.numeric(args[5])
wt_h=10
wt_j=10
gap_e=0.1
# gap_e=as.numeric(args[6])
# gap_o=as.numeric(args[7])

gap_o_list=c(-0.01, -0.1, -1,-5, -10,-20, -50, -100)

# seqfiles=list.files(fileseqs,pattern = "raw")

seqs=seqinr::read.fasta(fileseq)
mymrf=RNAmrf::read_mrf_renum(filemrf)


pbapply::pboptions(type="txt")
pbapply::pbsapply(gap_o_list, function(gap_o){
	a2b_1 = align_seq2mrf(seqs[[1]], mymrf, wt_h=wt_h,wt_j=wt_j,gap_ext=gap_e,gap_open=gap_o, debug = FALSE)
	a2b_2 = align_seq2mrf(seqs[[2]], mymrf, wt_j=wt_j,wt_h=wt_h,gap_ext=gap_e,gap_open=gap_o, debug = FALSE)
	seqaln = pair_a2b2aln(a2b_1, a2b_2, seqs = seqs)

	seqinr::write.fasta(seqaln,
	                  names = names(seqs),
	                  file.out = file.path(fileout, paste0(basename(fileseq),gap_o, ".mrfaln.fa")))


})

  

