#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(RNAmrf))
# suppressPackageStartupMessages(library(tidyr))
# suppressPackageStartupMessages(library(dplyr))




### funcs

mrfaln_seqs = function(seqs, mrf, init_method=1) {
  pbapply::pboptions(type="txt")
  aln_a2m=pbapply::pblapply(seqs,function(aseq){
    aseq.low=tolower(aseq)
    seq.enc= RNAmrf:::encode_seq(aseq.low)
    seq.int=seq.enc$seq_int_ungapped
    
    a2b=RNAmrf:::align_seq2mrf(aseq.low,
                               mrf = mrf,
                               gap_open = -3,
                               debug = FALSE,
                               iteration = 20,
                               init_method=init_method)
    # a2b_1b=a2b+1
    
    a2m=RNAmrf:::a2b2a2m(a2b,seq.int,mrflen = mrf$len)
    return(a2m)
  })
  
  return(aln_a2m)
}

### main

filemrf=args[1]
fileseqs=args[2]
fileout=args[3]

seqs=seqinr::read.fasta(fileseqs,forceDNAtolower=FALSE)

mymrf=read_mrf(filemrf)


out_a2m =mrfaln_seqs(seqs,mymrf)

seqinr::write.fasta(out_a2m,names=names(out_a2m), 
                    file.out=fileout)




