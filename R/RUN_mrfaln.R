#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(RNAmrf))

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))



read_mrf_renum = function(filemrf) {
  
  myalphabet = c("a", "u", "c", "g", "-")
  
  v1 = data.table::fread(cmd = paste("grep '^V'", filemrf))
  names(v1)[-1] = myalphabet
  
  w1 = data.table::fread(cmd = paste("grep  '^W'", filemrf))
  
  
  len=nrow(v1)
  len_a=length(myalphabet)
  
  # renumber MRF
  
  v1$i_ori=as.integer(gsub(".*\\[(.*?)\\].*","\\1",v1$V1))
  v1$i=1:nrow(v1)
  
  
  w1$i_ori=as.integer(gsub(".*\\[(.*?)\\]\\[(.*?)\\].*","\\1",w1$V1))
  w1$i=match(w1$i_ori,v1$i_ori)
  
  w1$j_ori=as.integer(gsub(".*\\[(.*?)\\]\\[(.*?)\\].*","\\2",w1$V1))
  w1$j=match(w1$j_ori, v1$i_ori)
  
  array_j = array(0, dim = c(len, len, len_a, len_a))
  
  for (m in 1:nrow(w1)) {
    id_i = w1$i[m]
    id_j = w1$j[m]
    
    mat = matrix(as.matrix(w1[m, 2:26]), 5, 5, byrow = TRUE)
    array_j[id_i, id_j, , ] = mat
    
  }
  
  ids = expand.grid(1:nrow(v1), 1:nrow(v1)) %>% filter(Var2 > Var1) %>% arrange(Var1)
  w1_mat = as.matrix(w1[, 2:(1 + length(myalphabet) * length(myalphabet))])
  
  mat_score = sapply(1:nrow(w1_mat), function(i) {
    tmpmat = matrix(w1_mat[i,], 5, 5)
    score = sqrt(sum(tmpmat[-5,-5] ^ 2))
    return(score)
  })
  ids$score = mat_score
  
  mat_mrf = matrix(0, nrow(v1), nrow(v1))
  
  mat_mrf[as.matrix(ids[, c(1, 2)])] = ids[, 3]
  mat_mrf[as.matrix(ids[, c(2, 1)])] = ids[, 3]
  
  mat_apc = RNAmrf:::APC_correction(mat_mrf)
  
  mrf = list(
    len = len,
    h = v1,
    j = w1,
    mat_mrf = mat_mrf,
    mat_apc = mat_apc,
    array_j=array_j
  )
  mrf_mat=RNAmrf:::mrf2mrf_mat(mrf)
  mrfh = as.matrix((mrf$h[, 2:6]))
  
  mrf$mrf_mat=mrf_mat
  mrf$mrf_h=mrfh
  
  return(mrf)
}

pair_a2b2aln = function(a2b_1, a2b_2, seqs) {
  
  last_idx = -1
  for (i in 1:length(a2b_1)) {
    if (a2b_1[i] == -1) {
      a2b_1[a2b_1 > last_idx] = a2b_1[a2b_1 > last_idx] + 1
      a2b_2[a2b_2 > last_idx] = a2b_2[a2b_2 > last_idx] + 1
      a2b_1[i] = last_idx + 1
    }
    last_idx = a2b_1[i]
  }
  
  last_idx = -1
  for (i in 1:length(a2b_2)) {
    if (a2b_2[i] == -1) {
      a2b_1[a2b_1 > last_idx] = a2b_1[a2b_1 > last_idx] + 1
      a2b_2[a2b_2 > last_idx] = a2b_2[a2b_2 > last_idx] + 1
      a2b_2[i] = last_idx + 1
    }
    last_idx = a2b_2[i]
  }
  
  a2b_1_1b=a2b_1+1
  a2b_2_1b=a2b_2+1
  seq_aln1 = character(max(a2b_1_1b, a2b_2_1b))
  seq_aln1[] = "-"
  seq_aln2 = character(max(a2b_1_1b, a2b_2_1b))
  seq_aln2[] = "-"
  
  seq_aln1[a2b_1_1b] = seqs[[1]]
  seq_aln2[a2b_2_1b] = seqs[[2]]
  
  return(list(
    seq_aln1=seq_aln1,
    seq_aln2=seq_aln2
  ))
}


### main

# filemrf="/home/tc/GIT/sandbox/mrftmp/mrf/5_8S_rRNA/tmp.sto.afa.mrf"
# fileout="/home/tc/GIT/sandbox/mrftmp/mrf/5_8S_rRNA/"
# fileseqs="/home/tc/GIT/sandbox/mrftmp/5_8S_rRNA/"

filemrf=args[1]
fileout=args[2]
fileseqs=args[3]


seqfiles=list.files(fileseqs,pattern = "raw")


mymrf=read_mrf_renum(filemrf)


pbapply::pboptions(type="txt")
dummy = pbapply::pblapply(seqfiles, function(seqfile) {
  seqs = seqinr::read.fasta(file.path(fileseqs, seqfile))
  
  a2b_1 = align_seq2mrf(seqs[[1]], mymrf, wt_h = 0.5, debug = FALSE)
  a2b_2 = align_seq2mrf(seqs[[2]], mymrf, wt_h = 0.5, debug = FALSE)
  seqaln = pair_a2b2aln(a2b_1, a2b_2, seqs = seqs)
  
  seqinr::write.fasta(seqaln,
                      names = names(seqs),
                      file.out = file.path(fileout, paste0(seqfile, ".mrfaln.fa")))
})
