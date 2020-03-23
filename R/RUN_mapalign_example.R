#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(RNAmrf)

mrf=read_mrf("/home/tc/GIT/RNA_PKALN/data/set1_rfa_cln/3VRS_A_RF01734/tmp.mrf")

ct_ref=RNASSP::read_ct("/home/tc/GIT/RNA_PKALN/data/set1_rfa_cln/3VRS_A_RF01734/3VRS_A.ct")

ct_onlypk=ct_ref

ct_onlypk$j=0

ct_onlypk$j[c(2:5, 14:17)]=ct_ref$j[c(2:5, 14:17)]

seqs=seqinr::read.fasta("/home/tc/GIT/RNA_PKALN/data/set1_rfa_cln/3VRS_A_RF01734/tmp.ungapped.msa")


pbapply::pboptions(type="txt")
a2b.l= pbapply::pbsapply(seqs[2:100], function(seq){
  a2b = align_seq2mrf(seq,mrf,wt_h = 0.5,debug = FALSE)
  # bench_ref= bench_a2b(a2b, seq,mrf, seq_ref=seqs[[1]], ct_ref=ct_ref)
  # bench=bench_a2b(a2b,seq,mrf,seq_ref = seqs[[1]],ct_ref = ct_onlypk)
  # bench["pairid_all"]=bench_ref["pairid"]
  # return(bench)
  return(a2b)
})

# save(a2b.l, file="a2b.l.rda")
writeLines(a2b.l, "./a2b.l.txt")

