#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


seqs = seqinr::read.fasta(args[1], forceDNAtolower = FALSE)

seqs.rf= lapply(seqs,function(seq){
  return(seq[seq %in% c("A","U","C","G","N","R","Y","M","-")])
})

seqmat=do.call(rbind,seqs.rf)

apsi=RNAmrf:::calc_apsi(seqmat)

# writeLines(apsi)
writeLines(paste0(basename(args[1])," apsi: ",format(apsi,digits = 3)))