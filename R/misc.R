read_pkaln_dir = function(testdir, debug=FALSE) {
  caseid = basename(testdir)
  
  pkaln=list()
  
  dirlist = list.dirs(testdir, recursive = FALSE)
  dirlist = dirlist[!grepl("input", dirlist)]
  
  if (debug){
    print(dirlist)
  }
  ## reference
  seqs.ref = seqinr::read.fasta(file.path(testdir, "input", paste0(caseid, ".a2m")), forceDNAtolower =
                                  FALSE)
  pkaln$seqidx.ref = RNAmrf:::msa_a2m2seqidx_all(seqs.ref)
  inputsto = paste0(testdir,"/",caseid,".sto")
  seqidx.aln.list = lapply(dirlist, function(dir) {
    bn = basename(dir)
    seqs.aln = seqinr::read.fasta(file.path(dir, paste0(caseid, ".", bn, ".a2m")), forceDNAtolower =
                                    FALSE)
    seqidx = RNAmrf:::msa_a2m2seqidx_all(seqs.aln)
    
    mrffile = file.path(dir, paste0(caseid, ".", bn, ".mrf"))
    
    
    ## convert mrf reference column to reference in the MSA
    if (file.exists(mrffile)){
      dfref=RNAmrf:::read_dfref(inputsto,mrffile)
      tmpref=matrix(0,nrow=nrow(seqidx),ncol=nrow(dfref))
      tmpref[,dfref$id_ref[dfref$id_mrf>0 & dfref$id_ref >0]]=
      seqidx[,dfref$id_mrf[dfref$id_mrf>0 & dfref$id_ref >0]]
      seqidx = tmpref
    }
    
    return(seqidx)
    
  })
  
  dfref_simple=RNAmrf:::read_dfref(inputsto)
  names(seqidx.aln.list) = basename(dirlist)

  pkaln$seqidx.aln.list=seqidx.aln.list
  
  pkaln$caseid=caseid
  pkaln$dfref=dfref_simple
  return(pkaln)
}

bench_pkaln = function(pkaln) {
  rslt.l = lapply((pkaln$seqidx.aln.list), function(seqidx.test) {
    RNAmrf:::bench_dfref(seqidx.test, pkaln$seqidx.ref, pkaln$dfref)
  })
  
  rslt.df = as.data.frame(do.call(rbind, rslt.l))
  
  rslt.df$method = names(pkaln$seqidx.aln.list)
  rslt.df$caseid=pkaln$caseid
  rslt.df
}

bench_pkaln_2 = function(pkaln) {
  rslt.l = lapply((pkaln$seqidx.aln.list), function(seqidx.test) {
    bench_dfref_2(seqidx.test, pkaln$seqidx.ref, pkaln$dfref)
  })
  
  rslt.df = as.data.frame(do.call(rbind, rslt.l))
  
  rslt.df$method = names(pkaln$seqidx.aln.list)
  rslt.df$caseid=pkaln$caseid
  rslt.df
}

## overall
bench_pkaln_all = function(pkaln) {
  tmp = lapply((pkaln$seqidx.aln.list), function(seqidx.test) {
    bench_dfref_all(seqidx.test, pkaln$seqidx.ref, pkaln$dfref)
  })
  
  rslt.df = as.data.frame(do.call(rbind, tmp))
  
  rslt.df$method = names(pkaln$seqidx.aln.list)
  rslt.df$caseid=pkaln$caseid
  rslt.df
}

bench_pkaln_2 = function(pkaln) {
  rslt.l = lapply((pkaln$seqidx.aln.list), function(seqidx.test) {
    bench_dfref_2(seqidx.test, pkaln$seqidx.ref, pkaln$dfref)
  })
  
  rslt.df = as.data.frame(do.call(rbind, rslt.l))
  
  rslt.df$method = names(pkaln$seqidx.aln.list)
  rslt.df$caseid=pkaln$caseid
  rslt.df
}



##

bench_pair_spfm = function(seqidx.test,
                           seqidx.ref,
                           pair, 
                           debug = FALSE,
                           raw = FALSE) {
  idx_A = which(pair > 0)
  idx_B = pair[which(pair > 0)]
  
  TPFN = sum(seqidx.ref[, idx_A] > 0)
  
  TPFP = sum(seqidx.test[, idx_A] > 0)
  
  TP = sum(seqidx.test[, idx_A] == seqidx.ref[, idx_A] &
             seqidx.test[, idx_B] == seqidx.ref[, idx_B] &
             seqidx.ref[, idx_A] > 0)
  
  sens = TP / TPFN
  prec = TP / TPFP
  mcc = sqrt(sens * prec)
  f1 = 2 * sens * prec / (sens + prec)
  
  # message(("FP | TP TPFN TPFP"))
  if (debug) {
    message(sprintf("%d |%d %d %d", (TPFP - TP), TP, TPFN, TPFP))
    
  }
  rslts = c(
    sens = sens,
    prec = prec,
    mcc = mcc,
    f1 = f1
  )
  rslts_raw = c(
    TP=TP,
    TPFP=TPFP,
    TPFN=TPFN
  )
  if (raw){
    return(rslts_raw)
  }
  else{
    return(rslts)
  }

}


bench_by_range_spfm = function(seqidx.test,
                               seqidx.ref,
                               range = NULL,
                               debug = FALSE,
                               raw = FALSE) {
  if (is.null(range))
    range = 1:ncol(seqidx.ref)
  
  TP = sum(seqidx.test[, range] == seqidx.ref[, range] &
             seqidx.ref[, range] > 0)
  TPFN = sum(seqidx.ref[, range] > 0)
  TPFP = sum(seqidx.test[, range] > 0)
  
  sens = TP / TPFN
  prec = TP / TPFP
  mcc = sqrt(sens * prec)
  f1 = 2 * sens * prec / (sens + prec)
  
  # message(("FP | TP TPFN TPFP"))
  if (debug) {
    message(sprintf("%d |%d %d %d", (TPFP - TP), TP, TPFN, TPFP))
    
  }
  rslts = c(
    sens = sens,
    prec = prec,
    mcc = mcc,
    f1 = f1
  )
  rslts_raw = c(
    TP=TP,
    TPFP=TPFP,
    TPFN=TPFN
  )
  if (raw){
    return(rslts_raw)
  }
  else{
    return(rslts)
  }
}

bench_dfref_2= function(seqidx.test, seqidx.ref, dfref){
  idx_select = RNAmrf:::get_idx_select(dfref)
  
  rslt=list()
  
  pair_all=dfref$id_ref_pair
  pair_pk=pair_all
  pair_pk[!dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d")]=0
  
  pair_nonpk=pair_all
  pair_nonpk[dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d")]=0
  
  
  rslt$col_all=bench_by_range_spfm(seqidx.test,seqidx.ref, idx_select$col_all)
  rslt$col_loop=bench_by_range_spfm(seqidx.test,seqidx.ref,idx_select$col_loop)
  rslt$pair_all=bench_pair_spfm(seqidx.test,seqidx.ref, pair_all)
  rslt$pair_pk=bench_pair_spfm(seqidx.test,seqidx.ref, pair_pk)
  rslt$pair_nonpk=bench_pair_spfm(seqidx.test,seqidx.ref, pair_nonpk)
  
  return(unlist(rslt))
}

bench_dfref_all = function(seqidx.test, seqidx.ref, dfref){
  idx_select = RNAmrf:::get_idx_select(dfref)
  
  rslt=list()
  
  pair_all=dfref$id_ref_pair
  pair_pk=pair_all
  pair_pk[!dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d")]=0
  
  pair_nonpk=pair_all
  pair_nonpk[dfref$ss %in% c("A", "a", "B", "b", "C", "c", "D", "d")]=0
  
  
  rslt$col_all=bench_by_range_spfm(seqidx.test,seqidx.ref, idx_select$col_all,raw=TRUE)
  rslt$col_loop=bench_by_range_spfm(seqidx.test,seqidx.ref,idx_select$col_loop,raw=TRUE)
  rslt$pair_all=bench_pair_spfm(seqidx.test,seqidx.ref, pair_all,raw=TRUE)
  rslt$pair_pk=bench_pair_spfm(seqidx.test,seqidx.ref, pair_pk,raw=TRUE)
  rslt$pair_nonpk=bench_pair_spfm(seqidx.test,seqidx.ref, pair_nonpk,raw=TRUE)
  
  return(unlist(rslt))
}

