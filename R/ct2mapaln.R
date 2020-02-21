#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# source("~/GIT/rna_ml/R/IO_funs.R")

read_ct=function(file){
  # ct=data.table::fread(file, skip=1,header=FALSE)
  ct = data.table::fread(cmd= sprintf("cat %s|sed 's/\t\t\t\t/\t/g'",file) ,header=FALSE)

  names(ct)=c(
    "i",
    "nt",
    "i-1",
    "i+1",
    "j",
    "ref"
  )

  return(ct)
}



if (length(args) <1){
  stop("Usage: ct2mapaln.R file_ct")
}



ct=read_ct(args[1])

len=nrow(ct)

l1=sprintf("LEN\t%d",len)

ct.pairs=ct[ct$j>0]

cons=paste0("CON\t",ct.pairs$i-1,"\t",ct.pairs$j-1,"\t",1,collapse="\n")

mapstr=paste0(l1,"\n",cons)

writeLines(mapstr)



