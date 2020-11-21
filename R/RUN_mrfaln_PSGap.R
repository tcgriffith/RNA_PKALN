#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(RNAmrf))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

### main

# filemrf="/home/tc/GIT/sandbox/mrftmp/mrf/5_8S_rRNA/tmp.sto.afa.mrf"
# fileout="/home/tc/GIT/sandbox/mrftmp/mrf/5_8S_rRNA/"
# fileseqs="/home/tc/GIT/sandbox/mrftmp/5_8S_rRNA/"

filemrf=args[1]
fileout=args[2]
fileseqs=args[3]

filemsa=args[4] # for gap profile, afa format


#######################


msa=seqinr::read.fasta(filemsa)
msa.mat=do.call(rbind,msa)
gap_profile=
  apply(msa.mat,2,function(col){
    log((1+sum(col=="-"))/length(msa))
  })

########################



seqfiles=list.files(fileseqs,pattern = "raw")


mymrf=read_mrf_renum(filemrf)


pbapply::pboptions(type="txt")
dummy = pbapply::pblapply(seqfiles, function(seqfile) {
  seqs = seqinr::read.fasta(file.path(fileseqs, seqfile))
  
  a2b_1=RNAmrf::align_seq2mrf_PSgap(
    seq = seqs[[1]],
    mrf = mymrf, 
    gap_open = gap_profile,
    debug=FALSE
  )
  
  a2b_2=RNAmrf::align_seq2mrf_PSgap(
    seq = seqs[[2]],
    mrf = mymrf, 
    gap_open = gap_profile,
    debug=FALSE
  )
  
  # a2b_1 = RNAmrf::align_seq2mrf(seqs[[1]], mymrf, wt_h = 1.0, debug = FALSE)
  # a2b_2 = align_seq2mrf(seqs[[2]], mymrf, wt_h = 1.0, debug = FALSE)
  seqaln = pair_a2b2aln(a2b_1, a2b_2, seqs = seqs)
  
  seqinr::write.fasta(seqaln,
                      names = names(seqs),
                      file.out = file.path(fileout, paste0(seqfile, ".mrfaln.fa")))
})
