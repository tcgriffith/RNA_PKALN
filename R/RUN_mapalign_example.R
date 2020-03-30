#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(RNAmrf)

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))



mrf=read_mrf("/home/tc/GIT/RNA_PKALN/data/set1_rfa_cln/3VRS_A_RF01734/tmp.mrf")

ct_ref=RNASSP::read_ct("/home/tc/GIT/RNA_PKALN/data/set1_rfa_cln/3VRS_A_RF01734/3VRS_A.ct")


seqs=seqinr::read.fasta("/home/tc/GIT/RNA_PKALN/data/set1_rfa_cln/3VRS_A_RF01734/tmp.ungapped.msa")

myseqs=seqs[2:10]


pbapply::pboptions(type="txt")
bench_100= pbapply::pblapply(myseqs, function(seq){
  a2b = align_seq2mrf(seq,mrf,wt_h = 0.5,debug = FALSE)
  # bench_ref= bench_a2b(a2b, seq,mrf, seq_ref=seqs[[1]], ct_ref=ct_ref)
  bench=bench_a2b(a2b,seq,mrf,seq_ref = seqs[[1]],ct_ref = ct_ref)
  # bench["pairid_all"]=bench_ref["pairid"]
  return(bench)
  # return(a2b)
})


bench_100_ref = pbapply::pblapply(myseqs, function(seq){
  exp_seq = encode_seq(seq)

  bench= bench_a2b(exp_seq$seq_int_ref-1, seq,mrf, seq_ref=seqs[[1]], ct_ref=ct_ref)
  return(bench)
})

bench_100.mapalign= data.frame(do.call(rbind,bench_100))
bench_100.mapalign$type = "mapalign"

bench_100.ref = data.frame(do.call(rbind,bench_100_ref))
bench_100.ref$type="ref"

bench_100.all = rbind(bench_100.mapalign,bench_100.ref)


bench_100.all$id = gsub("\\|.*","",rownames(bench_100.all))


bench_100.wider=
bench_100.all %>%
  pivot_wider(id_cols = id, names_from = type , values_from = c(1:5))

bench_100.wider %>%
  # arrange(pairid) %>%
  knitr::kable()
