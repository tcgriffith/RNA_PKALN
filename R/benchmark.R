bench_pairs= function(pred.pairs, ref.pairs){
  tp=sum(pred.pairs %in% ref.pairs)
  
  tpfp=length(pred.pairs)
  
  tpfn=length(ref.pairs)
  
  sens= tp/tpfn
  
  prec=tp/tpfp
  
  f1=2/(1/prec+1/sens)
  
  return(data.frame(
    prec=prec,
    sens=sens,
    f1=f1,
    stringsAsFactors=FALSE
  ))
}

bench_pairs_df=function(pred_pairs, ref_pairs,dim1){
  # pred_pairs = extract_aln(rmat)
  pred.pairs = id2_to_id1(pred_pairs$i, pred_pairs$j, dim = dim1)
  
  ref.pairs= id2_to_id1(ref_pairs$i,ref_pairs$j,dim=dim2)
  
  df=bench_pairs(pred.pairs, ref.pairs)
}


bench_seqid=function(seq,seqref){
  return(sum(seq==seqref)/length(seqref))
}


bench_pair=function(seq,seqref,ctref, debug=FALSE){
  
  npair=sum(ctref$j>0)
  
  pairs=paste0(seq[ctref$i[ctref$j>0]],seq[ctref$j[ctref$j>0]])
  
  pairs=toupper(pairs)
  
  if(debug){
    print(paste(pairs))
  }
  
  return(sum(pairs %in% RNASSP::energy2)/npair)
}

bench_aln=function(seq,seqref,ctref,debug=FALSE){
  seqid=bench_seqid(seq,seqref)
  pairid=bench_pair(seq,seqref,ctref,debug)
  return(c(
    seqid=seqid,
    pairid=pairid
  ))
}

bench_aln_list=function(seqlist,seqref,ctref,debug=FALSE){
  rslt.l=lapply(seqlist, function(seq){
    seqid = bench_seqid(seq, seqref)
    pairid=bench_pair(seq,seqref,ctref,debug)
    return(data.frame(seqid=seqid,pairid=pairid))
  })
  
  rslt.df=do.call(rbind, rslt.l)
  
  rslt.df$name=names(rslt.l)
  return(rslt.df)
  
}

