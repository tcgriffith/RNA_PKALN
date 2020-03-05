falign_R=function(score_mtx, rows,cols){
  max_sco = 0;
  # sco[rows+1,cols+1];
  gap_o=-1
  sco=matrix(0, rows+1,cols+1)
  for (i in 2:(rows+1)){
    for (j in 2:(cols+1)){
      A = sco[i-1,j-1] + score_mtx[i-1,j-1]
      D = sco[i-1,j]
      R = sco[i,j-1]
      
      if (A >= R) {
        if (A >= D) {
          sco[i, j] = A
        } else{
          sco[i,j] = D
        }
      }
      else{
        if (R >= D) {
          sco[i,j] = R
        } else{
          sco[i,j] = D
        }
      }
      
      if(sco[i,j] > max_sco){max_sco = sco[i,j]}
    }
  }
  return(max_sco)
  # return(sco)
}

align_R=function(score_mtx, gap_open=-1,gap_e=0.2,debug=FALSE){
  rows=nrow(score_mtx)
  cols=ncol(score_mtx)
  
  sco=matrix(0, rows+1,cols+1)
  label=matrix(0, rows+1,cols+1)
  max_sco=0
  a2b=integer(rows)
  a2b[]=-1
  max_i=0
  max_j=0
  
  for (i in 2:(rows+1)){
    
    for (j in 2:(cols+1)){
      A = sco[i-1,j-1] + score_mtx[i-1,j-1]
      D = sco[i-1,j]
      R = sco[i,j-1]
      if(label[i-1,j] == 1){D =D+ gap_open}else{D =D+ gap_open * gap_e}
      if(label[i,j-1] == 1){R =R+ gap_open}else{R =D+ gap_open * gap_e}
      # if(label[i-1,j] == 1){D =D+ gap_b[j-1]}else{D =D+ gap_b[j-1] * gap_e}
      # if(label[i,j-1] == 1){R =R+ gap_a[i-1]}else{R =D+ gap_a[i-1] * gap_e}
      
      if(A <= 0 && D <= 0 && R <= 0){label[i,j] = 0;sco[i,j] = 0;}
      else{
        if(A >= R){if(A >= D){label[i,j] = 1;sco[i,j] = A;}else{label[i,j] = 2;sco[i,j] = D;}}
        else{if(R >= D){label[i,j] = 3;sco[i,j] = R;}else{label[i,j] = 2;sco[i,j] = D;}}
        if(sco[i,j] > max_sco){max_i = i;max_j = j;max_sco = sco[i,j];}
      }
    }
  }
  
  i = max_i;
  j = max_j;
  
  while (1) {
    if (label[i,j] == 0) {
      break
    }
    else if (label[i,j] == 1) {
      a2b[i - 1] = j - 1
      i=i-1
      j=j-1
    }
    else if (label[i,j] == 2) {
      i=i-1
    }
    else if (label[i,j] == 3) {
      j=j-1
    }
  }
  
  if (debug){
    return(sco)
  }
  
  return(a2b)
  
}


read_mrf = function(filemrf) {
  myalphabet = c("a", "u", "c", "g", "-")
  len=nrow(v1)
  len_a=length(myalphabet)
  
  v1 = data.table::fread(cmd = paste("grep '^V'", filemrf))
  names(v1)[-1] = myalphabet
  
  w1 = data.table::fread(cmd = paste("grep  '^W'", filemrf))
  
  w1$i=as.integer(gsub(".*\\[(.*?)\\]\\[(.*?)\\].*","\\1",w1$V1))+1
  
  w1$j=as.integer(gsub(".*\\[(.*?)\\]\\[(.*?)\\].*","\\2",w1$V1))+1
  
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
  
  mat_apc = APC_correction(mat_mrf)
  mrf = list(
    len = len,
    h = v1,
    j = w1,
    mat_mrf = mat_mrf,
    mat_apc = mat_apc,
    array_j=array_j
  )
  
  return(mrf)
}
