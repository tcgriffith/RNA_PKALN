#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source(here::here("R/misc.R"))

if(length(args)<1){
  message("Usage: RUN pkalndir output.rds" )
  q("no",0)
}

dir=args[1]
fileout=args[2]


dirs=list.dirs(dir,recursive=FALSE,full.names=TRUE)

names(dirs)= basename(dirs)

pbapply::pboptions(type="txt")

pkaln.all = pbapply::pblapply(dirs,function(dir){
  try((read_pkaln_dir(dir)))
})

saveRDS(pkaln.all, file=fileout)

