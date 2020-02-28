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