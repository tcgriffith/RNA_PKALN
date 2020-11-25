read_pkaln_dir = function(testdir) {
  caseid = basename(testdir)
  
  pkaln=list()
  
  dirlist = list.dirs(testdir, recursive = FALSE)
  dirlist = dirlist[!grepl("input", dirlist)]
  
  seqs.ref = seqinr::read.fasta(file.path(testdir, "input", paste0(caseid, ".a2m")), forceDNAtolower =
                                  FALSE)
  pkaln$seqidx.ref = RNAmrf:::msa_a2m2seqidx_all(seqs.ref)
  
  seqidx.aln.list = lapply(dirlist, function(dir) {
    bn = basename(dir)
    seqs.aln = seqinr::read.fasta(file.path(dir, paste0(caseid, ".", bn, ".a2m")), forceDNAtolower =
                                    FALSE)
    RNAmrf:::msa_a2m2seqidx_all(seqs.aln)
  })
  
  dfref=RNAmrf:::read_dfref(paste0(testdir,"/",caseid,".sto"),
                                  file.path(testdir, "input", paste0(caseid, ".afa.mrf"))) 
  names(seqidx.aln.list) = basename(dirlist)
  ## fix rnamrf idx to idx_ref
  # if ("rnamrf" %in% names(seqidx.aln.list))
  {
    tmp = seqidx.aln.list$rnamrf
    # dfref=dfref
    tmpref=matrix(0,nrow=nrow(tmp),ncol=nrow(dfref))
    tmpref[,dfref$id_ref[dfref$id_mrf>0 & dfref$id_ref >0]]=tmp[,dfref$id_mrf[dfref$id_mrf>0 & dfref$id_ref >0]]
    seqidx.aln.list$rnamrf=tmpref
  }

  
  
  pkaln$seqidx.aln.list=seqidx.aln.list
  

  pkaln$caseid=caseid
  pkaln$dfref=dfref
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