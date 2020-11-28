#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


get_mrfss_from_sto_mrf=function(fsto, fmrf=NULL) {

  ss = RNAmrf::sto2ss(fsto) # seedss
  refs=RNAmrf:::sto2ref(fsto)
  refs.c=seqinr::s2c(refs)
  idx_ref = which(refs.c != ".")

  df = RNAmrf::ss2pairs(ss) # seed ss

  # idx_ref = which(df$ss != ".")
  df$id_ref = 0
  df$id_ref[idx_ref] = 1:length(idx_ref)


  df$id_ref_pair=df$id_ref[match(df$pair,df$id)]
  df$id_ref_pair[is.na(df$id_ref_pair)]=0
  df$id_mrf=df$id_ref

  if(!is.null(fmrf)){
    v1 = data.table::fread(cmd = paste("grep '^V'", fmrf))
    v1$i_ori=as.integer(gsub(".*\\[(.*?)\\].*","\\1",v1$V1))+1
    v1$i=1:nrow(v1)
    df$id_mrf=0
    df$id_mrf[v1$i_ori]=v1$i
    
    df$id_mrf_pair=df$id_mrf[match(df$pair,df$id)]
    df$id_mrf_pair[is.na(df$id_mrf_pair)]=0
  }
  dfref=df[df$id_mrf>0,]
  
  ss_mrf=dfref$ss
  
  # clean one sided pairs
  ss_mrf[!(dfref$pair %in% dfref$id | dfref$pair==0)] ="."

  # dfref = df[df$id_ref > 0, ]
  # dfref=df
  
  ss_mrf.s=seqinr::c2s(as.character(ss_mrf))
  return(ss_mrf.s)
  
}

###
fsto=args[1]
fmrf=args[2]
fout=args[3]

ssout=get_mrfss_from_sto_mrf(fsto,fmrf)

ssout.full=paste0("#=GC SS_cons ",ssout)
writeLines(ssout.full,fout)

