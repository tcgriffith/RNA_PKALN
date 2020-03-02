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