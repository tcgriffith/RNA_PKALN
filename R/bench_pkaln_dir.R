#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source(here::here("R/misc.R"))

if(length(args)<1){
  message("Usage: RUN pkalndir" )
  q("no",0)
}

# testdir="/home/tc/GIT/RNA_PKALN/data/synth/"

testdir=args[1]

mypkaln=read_pkaln_dir(testdir)

bench_pkaln(mypkaln)
