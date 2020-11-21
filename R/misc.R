ss2pairs = function(ss) {
  pairs = integer(length = length(ss))
  
  stack1 = integer() # <>
  
  stack2 = integer() # Aa
  
  stack3= integer() # ()
  
  for (i in 1:length(ss_c)) {
    if (ss[i] == "<") {
      stack1 = append(stack1, i)
    }
    if (ss[i] == ">") {
      pairs[tail(stack1, 1)] = i
      pairs[i] =tail(stack1, 1)
      stack1 = head(stack1, -1)
    }
    if (ss[i] == "A") {
      stack2 = append(stack2, i)
    }
    if (ss[i] == "a") {
      pairs[tail(stack2, 1)] = i
      pairs[i] =tail(stack2, 1)
      stack2 = head(stack2, -1)
    }
    if (ss[i] == "(") {
      stack3 = append(stack3, i)
    }
    if (ss[i] == ")") {
      pairs[tail(stack3, 1)] = i
      pairs[i] =tail(stack3, 1)
      stack3 = head(stack3, -1)
    }
  }
  return(pairs)
}

bench_alna2b = function(a2bs, a2brefs) {
  aln_sens = sapply(1:length(a2bs), function(i) {
    # a2bref = seqs[[i]]$seq_int_ref_renum
    a2bref=a2brefs[[i]]
    a2btest = a2bs[[i]]
    idx = which(a2btest > 0)
    ncols = sum(a2bref > 0)
    sens = sum(a2btest[idx] == a2bref[idx]) / ncols
    return(sens)
  })
  return((aln_sens))
  # return(aln_sens)
}


rng_seq=function(pairs=NULL,len=25){
  # myseq=character(25)
  myseq=sample(c("a","u","c","g","-"),25,replace = TRUE)
  
  
  ## add pairs
  i=which(pairs>0)
  j=pairs[which(pairs>0)]
  # message(i)
  
  rng_bp=tolower(sample(c("AU", "GC", "CG", "UA"),length(i),replace = TRUE))
  bp_1= gsub("^(.)(.)","\\1",rng_bp)
  bp_2= gsub("^(.)(.)","\\2",rng_bp)
  myseq[i]=bp_1
  myseq[j]=bp_2
  
  return(myseq)
}


a2b2a2m = function(a2b, seq,len_seq=25) {
  
  seq_enc = RNAmrf::encode_seq(seq)
  
  a2b_1 = a2b
  
  last_idx = 0
  for (i in 1:length(a2b)) {
    if (a2b_1[i] == 0) {
      a2b_1[a2b_1 > last_idx] = a2b_1[a2b_1 > last_idx] + 1
      a2b_1[i] = last_idx + 1 # fill unaligned
    }
    last_idx = a2b_1[i]
  }
  
  newseq = seq_enc$seq_ungapped
  newseq[a2b > 0] = toupper(newseq[a2b > 0])
  
  a2mlen = len_seq + sum(a2b == 0)
  
  # message(a2b_1)
  
  seq_a2m = character(a2mlen)
  seq_a2m[] = "-"
  seq_a2m[a2b_1] = newseq
  return(seq_a2m)
}

