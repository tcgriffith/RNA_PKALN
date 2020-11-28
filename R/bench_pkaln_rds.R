#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) <1){
  message("Usage: RUN pkaln.rds")
}

suppressPackageStartupMessages(library(dplyr))
options(warn=-1,dplyr.summarise.inform=FALSE)
source(here::here("R/misc.R"))


# testdir=args[1]

filepkaln=args[1]

pkaln.all=readRDS(filepkaln)

rslt.all = lapply(pkaln.all,function(pkaln){
  try(bench_pkaln(pkaln))
})

rslt.df.s = do.call(rbind,rslt.all)

print(rslt.df.s,digits = 2)

# print(rslt.all)

# mypkaln=read_pkaln_dir(testdir)

# bench_pkaln(mypkaln)


# filects=args[1]


# 
# cts.all.bind=ctsall_bind(cts.all)
# 
# rslt.df=bench_cts(cts.all.bind)
# 
# rslt.df.s = rslt.df %>% arrange(type) 

