#!/usr/bin/env Rscript

library(RNAmrf)
source(here::here("R/misc.R"))

pairs=rep(0,25)

pairs[1:4]=15:12
pairs[6:9]=22:19
# pairs

ss=rep(".",25)
ss[1:4]="("
ss[15:12]=")"
ss[6:9]="A"
ss[22:19]="a"

# seqinr::c2s(ss)
ss_s=seqinr::c2s(ss)

rng_seq=function(pairs=NULL,len=25){
  # myseq=character(25)
  myseq=sample(c("a","u","c","g","-"),25,replace = TRUE)
  
  
  ## add pairs
  i=which(pairs>0)
  j=pairs[which(pairs>0)]
  # message(i)
  
  rng_bp=tolower(sample(c("AU", "GC", "CG", "UA"),length(i),replace = TRUE))
  bp_1= gsub("^(.)(.)","\\1",rng_bp)
  bp_2= gsub("^(.)(.)","\\2",rng_bp)
  myseq[i]=bp_1
  myseq[j]=bp_2
  
  return(myseq)
}

seq1=rng_seq(pairs=pairs,len=25)

set.seed(42)
seqtmp=lapply(1:1000,function(i){
  return(rng_seq(pairs=pairs,len=25)) 
})

seqnames=paste0("seq",1:1000)

dir.create(here::here("data/synth2/"))

##ss
seqinr::write.fasta(ss_s,names="refSS",file.out=here::here("data/synth2/seq1000_SS.fasta"))

seqinr::write.fasta(seqtmp,names=seqnames,file.out=here::here("data/synth2/seq1000_gaps.fasta"))