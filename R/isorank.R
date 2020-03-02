gg_mat = function(mat) {
  # df=as.data.frame(which(mat>0,arr.ind=TRUE))
  
  df = as.data.frame(which(!is.na(mat) > 0, arr.ind = TRUE))
  
  df$val = mat[which(!is.na(mat) > 0)]
  # df=as.data.frame(as.table(mat))
  ggplot(df) + geom_tile(aes(x=col, y=row, fill = val)) + coord_fixed()
}

extract_aln = function(R, debug=FALSE) {
  Rtmp = R
  maplist=list()
  
  for (i in 1:nrow(Rtmp)) {
    max_R = max(Rtmp)
    if (max_R > 0.01) {
      
      
      
      idmax = which.max(Rtmp)
      ids = id1_to_id2(idmax, nrow(Rtmp))
      maplist=append(maplist,
                     values=list(
                       data.frame(
                         i=ids$i[[1]],
                         j=ids$j[[1]],
                         v=Rtmp[idmax]
                       )
                     ))
      if(debug){
        message(sprintf("%d %d %d %f", i, ids$i[[1]], ids$j[[1]], Rtmp[idmax]))
      }
      # Rtmp=Rtmp[-ids["i"],-ids["j"]]
      Rtmp[ids$i[[1]], ] = -1
      Rtmp[, ids$j[[1]]] = -1
      # gg=gg_mat(Rtmp)
      # plot(gg)
    }
    
  }
  
  return(do.call(rbind,maplist))
  # return(maplist)
}

id1_to_id2=function(id,dim){
  # j = id %/% (dim +1)+1
  i = id %% (dim)
  j=ceiling(id/dim)
  i[i==0]=dim
  return(data.frame(i=i,j=j))
}

id2_to_id1=function(i,j,dim){
  return((j-1)*dim+i)
}



# not suitable for large matrices
# matA2matR=function(A,nrow){
#   A.eig=eigen(A)
#   R = matrix(Re(A.eig$vectors[,1]),nrow=nrow)
#   return(R)
# }

Rv2Rmat=function(Rv,nrow,ncol){
  Rmat=matrix(Rv, nrow,ncol)
}

graph2mat=function(graph){
  mat=matrix(0,nrow=nrow(graph$nodes),ncol=nrow(graph$nodes))
  
  mat[as.matrix(graph$edges[,c("from","to")])]=1
  mat[as.matrix(graph$edges[,c("to","from")])]=1
  return(mat)
}

mat2graph=function(mat,id_renum=0){
  
  nodes=data.frame(
    id=1:nrow(mat)+id_renum,
    label=row.names(mat),
    stringsAsFactors=FALSE
  )
  
  edges.mat=which(mat>0,arr.ind=TRUE)
  
  edges=as.data.frame(edges.mat+id_renum)
  colnames(edges)= c("from","to")
  
  edges=edges[edges$from < edges$to,]
  
  mygraph=list()
  mygraph$nodes=nodes
  mygraph$edges=edges
  return(mygraph)
}

id2_to_id1=function(i,j,dim){
  return((j-1)*dim+i)
}

get_Rmatv = function(A, l1, l2,iteration=1000,debug=FALSE) {
  Rv = matrix(0, nrow = l1 * l2, ncol = iteration)
  
  R=runif(l1*l2)
  
  Rv[, 1] = R
  
  for (i in 2:iteration) {
    Rvnew = A %*% Rv[, i - 1]
    Rvnew = as.matrix(Rvnew)
    Rv[, i] = Rvnew / norm(Rvnew, "F")
    
    lastdist=dist(rbind(Rv[, i], Rv[, i - 1]))
    


    if(debug){
      if (i %% 100 ==0) message(i, " ", lastdist)
    }

  }
  if(debug){
    message("##  ", i, " ", lastdist)
  }
  
  
  # Rmat = Rv2Rmat(, nrow = l1, ncol = l2)
  
  return(Rv)
  
  
}

get_matA=function(mat1,mat2){
  l1=nrow(mat1)
  l2=nrow(mat2)
  
  edges1 = which(mat1 > 0)
  edges2 = which(mat2 > 0)
  
  edgepairs = expand.grid(edges1, edges2)
  
  u_contacts = colSums(mat1)
  v_contacts = colSums(mat2)
  
  
  edgepairs1 = arrayInd(edgepairs$Var1, .dim = c(l1, l1))
  colnames(edgepairs1) = c("i", "u")
  edgepairs2 = arrayInd(edgepairs$Var2, .dim = c(l2, l2))
  colnames(edgepairs2) = c("j", "v")
  
  edgepairsall = cbind(edgepairs1, edgepairs2)
  
  edgepairsall.df =
    edgepairsall %>%
    as.data.frame() %>%
    mutate(
      id_row = id2_to_id1(i, j, l1),
      id_col = id2_to_id1(u, v, l1),
      V = 1 / (u_contacts[u] * v_contacts[v])
    )
  
  A_sparse = Matrix::sparseMatrix(
    i = edgepairsall.df$id_row,
    j = edgepairsall.df$id_col,
    x = edgepairsall.df$V,
    dims = c(l1 * l2, l1 * l2)
  )
  return(A_sparse)
}

run_isorank = function(mat1, mat2,iteration=1000,debug=FALSE) {
  
  l1=nrow(mat1)
  l2=nrow(mat2)
  
  A_sparse=get_matA(mat1,mat2)
  
  message("## A_sparse constructed, calculating R")
  
  Rmatv=get_Rmatv(A_sparse,l1,l2,iteration,debug)
  message("## Rmat extracted, finished")  
  
  Rmat=matrix(Rmatv[, iteration], nrow=l1,ncol=l2)
  # return(Rmat)
  return(Rmat)
}


get_A_nb = function(mat1, mat2) {
  l1 = nrow(mat1)
  l2 = nrow(mat2)
  
  # A_nb = Matrix::sparseMatrix(i=NULL,j=NULL,dims=c(l1*l2,l1*l2))
  
  id_col=integer()
  id_row=integer()
  
  for (i in 1:l1) {
    for (j in 1:l2) {
      idx_ij = id2_to_id1(i, j, l1)
      idx_imjm = id2_to_id1(i - 1, j - 1, l1)
      idx_ipjp = id2_to_id1(i + 1, j + 1, l1)
      
      if (idx_imjm > 0) {
        # A_nb[idx_ij, idx_imjm] = 1 / 2
        id_row=append(id_row, idx_ij)
        id_col=append(id_col,idx_imjm)
      }
      if (idx_ipjp < l1 * l2) {
        id_row=append(id_row,idx_ij)
        id_col=append(id_col,idx_ipjp)
        # A_nb[idx_ij, idx_ipjp] = 1 / 2
      }
    }
  }
  
  A_nb = Matrix::sparseMatrix(i=id_row,j=id_col,x=1/2,dims=c(l1*l2,l1*l2))
  
  return(A_nb)
}

get_Aseq = function(mat1, mat2) {
  seq1=colnames(mat1)
  seq2=colnames(mat2)
  
  
  matB = matrix(0, nrow = length(seq1), ncol = length(seq2))
  
  for (i in 1:nrow(matB)) {
    for (j in 1:ncol(matB)) {
      if (seq1[i] == seq2[j]) {
        matB[i, j] = 1
      }
    }
  }
  
  matB = matB / sqrt(sum(matB ^ 2))
  
  vB = as.vector(matB)
  
  # vB=vB/sum(vB)
  
  A_seq = vB %*% t(vB)
  
  A_seq.sparse=Matrix::Matrix(A_seq,sparse=TRUE)
  
  return(A_seq.sparse)
}



ct2mat = function(ct, add_backbone=FALSE) {
  l1 = nrow(ct)
  mat1 = matrix(0, l1, l1)
  
  rownames(mat1) = ct$i
  
  colnames(mat1) = ct$nt
  
  linkpair = ct[ct$j > 0, c("i", "j")]
  names(linkpair) = c("from", "to")
  
  linknb = ct[ct$`i+1` > 0, c("i", "i+1")]
  names(linknb) = c("from", "to")
  
  mat1[as.matrix(linkpair[, c(1, 2)])] = 1
  mat1[as.matrix(linkpair[, c(2, 1)])] = 1
  
  if(add_backbone){
    mat1[as.matrix(linknb[, c(1, 2)])] = 1
    mat1[as.matrix(linknb[, c(2, 1)])] = 1
  }
  
  # mat1[as.matrix(linknb[, c(1, 2)])] = 1
  # mat1[as.matrix(links[, c(2, 1)])] = 1
  return(mat1)
}


run_isorank_bpnbseq = function(mat1,
                               mat2,
                               alpha = 0.3,
                               beta = 0.3,
                               iteration = 1000,
                               debug = FALSE) {
  # alpha = alphh
  # beta = .4
  l1 = nrow(mat1)
  l2 = nrow(mat2)
  
  A_bp = get_matA(mat1, mat2)
  A_nb = get_A_nb(mat1, mat2)
  A_seq = get_Aseq(mat1, mat2)
  
  
  A_all = alpha * A_bp + beta * A_nb + (1 - alpha - beta) * A_seq
  # A_all = alpha * A_bp + (1 - alpha) * A_nb
  # rmat = do_it(A_all, l1, l2, iteration = 10000, debug = FALSE)
  
  Rmatv=get_Rmatv(A_all,l1,l2,iteration,debug)
  Rmat=matrix(Rmatv[, iteration], nrow=l1,ncol=l2)
  # gg_mat(rmat)
  return(Rmat)
}

run_isorank_A3 = function(A_bp,
                          A_nb,
                          A_seq,
                          alpha = 0.3,
                          beta = 0.3,
                          iteration = 1000,
                          debug = FALSE) {
  # alpha = alphh
  # beta = .4
  l1 = nrow(mat1)
  l2 = nrow(mat2)
  
  # A_bp = get_matA(mat1, mat2)
  # A_nb = get_A_nb(mat1, mat2)
  # A_seq = get_Aseq(mat1, mat2)
  
  
  A_all = alpha * A_bp + beta * A_nb + (1 - alpha - beta) * A_seq
  # A_all = alpha * A_bp + (1 - alpha) * A_nb
  # rmat = do_it(A_all, l1, l2, iteration = 10000, debug = FALSE)
  
  Rmatv = get_Rmatv(A_all, l1, l2, iteration, debug)
  Rmat = matrix(Rmatv[, iteration], nrow = l1, ncol = l2)
  # gg_mat(rmat)
  return(Rmat)
}



run_isorank_bpnb = function(mat1,
                            mat2,
                            alpha = 0.5,
                            iteration = 1000,
                            debug = FALSE) {
  l1=nrow(mat1)
  l2=nrow(mat2)
  A_bp = get_matA(mat1, mat2)
  A_nb = get_A_nb(mat1, mat2)
  
  A_all = alpha * A_bp + (1 - alpha) * A_nb
  
  Rmatv=get_Rmatv(A_all,l1,l2,iteration,debug)
  Rmat=matrix(Rmatv[, iteration], nrow=l1,ncol=l2)
  # gg_mat(rmat)
  return(Rmat)
  
}


