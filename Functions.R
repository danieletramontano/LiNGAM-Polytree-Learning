##LIBRARIES
library(igraph)
library(Matrix)
library(CompareCausalNetworks)
library(kernlab)
library(bnlearn)
library(pcalg)


###Chow-Liu+Pairwise Orientation with K=4

#Takes as input the data matrix x, the covariance matrix S and the unoriented edge e=(i,j)
#Return 1 if i->j is the right orientation and 2 otherwise 
edgeorient_third_fourth_cum<-function(x,S,e){
  i<-e[1]
  j<-e[2]
  M1<-matrix(c(S[i,i],tmom(x,i,i,i),tmom(x,i,i,j),fmom(x,i,i,i,i)-3*S[i,i]^2,fmom(x,i,i,i,j)-3*S[i,i]*S[i,j],fmom(x,i,i,j,j)-S[i,i]*S[j,j]-2*S[i,j]^2,S[j,i],0,tmom(x,j,i,j),0,0,fmom(x,i,j,j,j)-3*S[j,j]*S[i,j]),nrow=2,byrow=TRUE)
  M1[2,2]<-M1[1,3]
  M1[2,4]<-M1[1,5]
  M1[2,5]<-M1[1,6]
  M2<-matrix(c(M1[2,1],M1[1,3],M1[2,3],M1[1,5],M1[1,6],M1[2,6],S[j,j],M1[2,3],tmom(x,j,j,j),M1[1,6],M1[2,6],fmom(x,j,j,j,j)-3*S[j,j]^2),nrow=2,byrow=TRUE)
  a1<-norm(as(c(det(M1[,c(1,2)]),det(M1[,c(1,3)]),det(M1[,c(1,4)]),det(M1[,c(1,5)]),det(M1[,c(1,6)])),"matrix"),type = "2")
  a2<-norm(as(c(det(M2[,c(1,2)]),det(M2[,c(1,3)]),det(M2[,c(1,4)]),det(M2[,c(1,5)]),det(M2[,c(1,6)])),"matrix"),type = "2")
  return(which.min(c(a1,a2)))
}


#Takes as input the data matrix x, the covariance matrix S and the list of unoriented edges E
#Returns adjacency matrix I.
Oneedgeordering_third_fourth_cum<-function(x,S,E){
  I<-Matrix(nrow=length(S[1,]),ncol=length(S[1,]),data=0,sparse=TRUE)
  for(i in c(1:length(E[,1]))){
    m<-edgeorient_third_fourth_cum(x,S,E[i,])
    if(m==1){
      I[E[i,1],E[i,2]]<-1
    }
    else{
      I[E[i,2],E[i,1]]<-1
    }
  }
  return(I)
}

#Takes as input the data matrix x.
#Return the adjacency matrix I.
Onedgepolytree_third_fourth_cum<-function(x){
  S<-cov(x)
  n<-dim(S)[1]
  C<-cov2cor(S)
  w<--abs(as.vector(C[t(upper.tri(C))]))
  g<-make_full_graph(n)
  m<-mst(g,w)
  E<-get.edgelist(m)
  I<-Oneedgeordering_third_fourth_cum(x,S,E)
  return(I)
}



##Chow-Liu + bivariateANM
#Takes as input the data matrix x, the covariance matrix S and the list of unoriented edges E.
#Returns adjacency matrix I.
Oneedgeordering_anm<-function(x,S,E){
  I<-Matrix(nrow=length(S[1,]),ncol=length(S[1,]),data=0,sparse=TRUE)
  for(i in c(1:length(E[,1]))){
    m<-getParents(x[,E[i,]],method = "bivariateANM",sparse=TRUE)
    if(m[1,2]==1){
      I[E[i,1],E[i,2]]<-1
    }
    else if(m[2,1]==1){
      I[E[i,2],E[i,1]]<-1
    }
  }
  return(I)
}

#Takes as input the data matrix x.
#Return the adjacency matrix I.
Onedgepolytree_anm<-function(x){
  S<-cov(x)
  n<-dim(S)[1]
  C<-cov2cor(S)
  w<--abs(as.vector(C[t(upper.tri(C))]))
  g<-make_full_graph(n)
  m<-mst(g,w)
  E<-get.edgelist(m)
  I<-Oneedgeordering_anm(x,S,E)
  return(I)
}

###Chow-Liu+PT with K=4

#Takes as input the data matrix x, and the threshold. 
#Returns the adjacency matrix I, after performing the algorithm PT.
triplet_third_fourth_cum<-function(x,thres){
  S<-cov(x)
  n<-dim(S)[1]
  C<-cov2cor(S)
  w<--abs(as.vector(C[t(upper.tri(C))]))
  g<-make_full_graph(n)
  m<-mst(g,w)
  E<-get.edgelist(m)
  V<-c(1:n)
  LIST<-triplets(E,C,V,thres)
  Ulist<-LIST$Ulist
  Olist<-LIST$Olist
  Olist1<-RecursiveOneEdge_third_fourth_cum(x,C,S,Ulist,Olist,numeric())
  I<-edgelist_toadjmatrix(Olist1)
  return(I)
}


#Takes as imput the list of unoriented edges E, the correlation matrix C, the list of vertices V, and a threshold thres. 
#Returns the list of oriented edges O, and the list of unoriented edges E.

triplets<-function(E,C,V,thres){
  m<-length(V)
  V_v<-V
  UE<-E
  O<-numeric()
  c<-matrix(0,nrow = m,ncol=1)
  vcheck<-numeric()
  repeat{
    if(length(V_v)==0){
      LIST<-list(Olist=O,Ulist=UE)
      return(LIST)
      break
    }
    update<-FALSE
    for(i in c(1:length(V_v))){
      vupdate<-FALSE
      vcheck<-numeric()
      v<-V_v[i]
      if(length(UE)>2){
        Lu1<-which(UE[,1]==v)
        Lu2<-which(UE[,2]==v)
      }
      else{
        Lu1<-which(UE[1]==v)
        Lu2<-which(UE[2]==v)
      }
      L<-append(Lu1,Lu2)
      l<-length(L)
      if(l==1){
        vcheck<-append(vcheck,i)
      }
      if(length(O)>0){
        Li<-which(O[,2]==v)
      }
      else{
        Li<-numeric()
      }
      if((length(Li)>0)&(l>0)){
        e<-O[Li[1],]
        O<-withincoming(C,e,v,Lu1,Lu2,O,UE,thres)
        vcheck<-append(vcheck,i)
        update<-TRUE
        vupdate<-TRUE
      }
      else{
        if(l>1){
          d<-onlyundirected(C,v,UE,Lu1,Lu2,thres)
          col<-which(d==1)
          if(length(col)>0){
            a<-col[1]%%l
            if(a==0){
              a<-l
            }
            d<-d+t(d)
            d<-d+diag(dim(d)[1])
            
            d<-d[a,]
            if(a<=length(Lu1)){
              e<-UE[Lu1[a],2]
              O<-withlist(Lu1,Lu2,O,UE,d,thres)
            }
            else{
              e<-UE[Lu2[a],1]
              O<-withlist(Lu1,Lu2,O,UE,d,thres)
            }
            vcheck<-append(vcheck,i)
            update<-TRUE
            vupdate<-TRUE
          }
        }
      }
      if(vupdate){
        if(length(UE)>2){
          UE<-UE[-L,]
        }
        else{
          LIST<-list(Olist=O,Ulist=numeric())
          return(LIST)
          break
        }
      }
    }
    if(update==FALSE){
      LIST<-list(Olist=O,Ulist=UE)
      return(LIST)
      break
    }
    V_v<-V_v[-vcheck]
  }
}



#Takes as input the correlation matrix C, a possible unshielded collider i-j-k, and the threshold.
#Returns 1 the triple forms a collider and 0 otherwise.
tripletcomp<-function(C,i,j,k,thres){
  if(abs(C[i,k])<thres){
    r<-1
  }
  else{
    r<-0
  }
  return(r)
}



#Takes as input the correlation matrix C, an oriented edge i->v, a vertex v, two lists of edges Lu1/Lu2 having v as a vertex, the lists of oriented and unoriented edges 
#and the threshold. 
#Returns the updated lists of oriented.
withincoming<-function(C,e,v,Lu1,Lu2,O,UE,thres){
  d<-matrix(0,nrow = length(Lu1)+length(Lu2),ncol=1) 
  if(length(Lu1)>0){
    for(j in c(1:length(Lu1))){
      if(length(UE)>2){
        d[j]<-tripletcomp(C,e[1],v,UE[Lu1[j],2],thres)
      }
      else{
        d[j]<-tripletcomp(C,e[1],v,UE[2],thres)
      }
    }
  }
  if(length(Lu2)>0){
    for(j in c(1:length(Lu2))){
      if(length(UE)>2){
        d[length(Lu1)+j]<-tripletcomp(S,e[1],v,UE[Lu2[j],1],thres)
      }
      else{
        d[length(Lu1)+j]<-tripletcomp(S,e[1],v,UE[1],thres)
      }
    }
  }
  O<-withlist(Lu1,Lu2,O,UE,d)
  return(O)
}

#Takes as input two correlation matrix C, a  vertex v, the two lists Lu1/Lu2 of edges having v as a vertex, and the threshold.
#Returns a matrix d, in which d[i,j]=1 if the edges in position i and j forms an unshielded collider and 0 otherwise.

onlyundirected<-function(C,v,UE,Lu1,Lu2,thres){
  L<-append(Lu1,Lu2)
  l<-length(L)
  d<-matrix(0,nrow =l,ncol=l)
  for(i in c(1:(l-1))){
    for(j in c((i+1):l)){
      if(i<=length(Lu1)){
        count_i<-2
      }
      else{
        count_i<-1
      }
      if(j<=length(Lu1)){
        count_j<-2
      }
      else{
        count_j<-1
      }
      d[i,j]<-tripletcomp(C,UE[L[i],count_i],v,UE[L[j],count_j],thres)
    }
  }
  return(d)
}



#Takes as input two lists Lu1/Lu2 sharing a given vertex, the lists of oriented and unoriented edges and the threshold, and the list of unshielded colliders
#with center the vertex. 
#Returns the updated list of oriented edges

withlist<-function(Lu1,Lu2,O,UE,d,thres){
  for(j in c(1:length(d))){
    if(j<=length(Lu1)){
      if(d[j]==0){
        if(length(UE)>2){
          O<-rbind(O,UE[Lu1[j],])  
        }
        else{
          O<-rbind(O,UE)  
        }
      }
      else{
        if(length(UE)>2){
          O<-rbind(O,rev(UE[Lu1[j],]))
        }
        else{
          O<-rbind(O,rev(UE))
        }
      }
    }
    else{
      if(d[j]==0){
        if(length(UE)>2){
          O<-rbind(O,rev(UE[Lu2[j-length(Lu1)],]))
        }
        else{
          O<-rbind(O,rev(UE))
        }
      }
      else{
        if(length(UE)>2){
          O<-rbind(O,UE[Lu2[j-length(Lu1)],])
        }
        else{
          O<-rbind(O,UE)
        }
      }
    }
  }
  return(O)
}

#Takes as input the data matrix x, the covariance matrix S and the lists of oriented and unoriented edges.
#Returns the list of Oriented edges
RecursiveOneEdge_third_fourth_cum<-function(x,S,E,O){
  TOUSE<-numeric()
  n<-dim(S)[1]
  repeat{
    if(length(O)!=(2*n-2)){
      if(length(TOUSE)==0){
        E<-matrix(E,nrow=length(E)/2)
        orient<-edgeorient_third_fourth_cum(x,S,E[1,])
        if(orient==1){
          o<-E[1,]
        }
        else{
          o<-rev(E[1,])
        }
        O<-rbind(O,o)
        E<-E[-1,]
        E<-matrix(E,nrow=length(E)/2)
        if(length(E)==0){
          return(O)
          break
        }
      }
      else{
        TOUSE<-matrix(TOUSE,nrow=length(TOUSE)/2)
        o<-TOUSE[1,]
        TOUSE<-TOUSE[-1,]
      }
      if(length(E)!=0){
        E<-matrix(E,nrow=length(E)/2)
        l1<-which(E[,1]==o[2])
        l2<-which(E[,2]==o[2])
        for(i in l2){
          E[i,]<-rev(E[i,])
        }
        l<-append(l1,l2)
        E_o<-matrix(E[l,],nrow=length(E[l,])/2)
        if(length(E_o)==0){
        }
        else{
          E<-E[-l,]
          for(i in c(1:(length(E_o)/2))){
            O<-rbind(O,E_o[i,])
            TOUSE<-rbind(TOUSE,E_o[i,])
          }
        }
      }
    }
    else{
      return(O)
      break
    }
  }
}




###Chow-Liu+TP with K=4

#Takes as input the data matrix x and a threshold.
#Returns the asjacency matrix I, after performing TP.
tp_third_fourth_cum<-function(x,thres){
  S<-cov(x)
  n<-dim(S)[1]
  C<-cov2cor(S)
  w<--abs(as.vector(C[t(upper.tri(C))]))
  g<-make_full_graph(n)
  m<-mst(g,w)
  E<-get.edgelist(m)
  O<-tmomcor_third_fourth_cum(x,C,S,E,thresh)
  I<-edgelist_toadjmatrix(O)
  return(I)
}


#Takes as input the data matrix x, the correlation matrix C, the covariance matrix S, the unoriented list of edges E and a threshold.
#Returns the oriented list of edges 0.
tmomcor_third_fourth_cum<-function(x,C,S,E,thres){
  TOUSE<-numeric()
  n<-dim(C)[1]
  O<-numeric()
  repeat{
    if(length(O)!=(2*n-2)){
      if(length(TOUSE)==0){
        E<-matrix(E,nrow=length(E)/2)
        orient<-edgeorient_third_fourth_cum(x,S,E[1,])
        if(orient==1){
          o<-E[1,]
        }
        else{
          o<-rev(E[1,])
        }
        O<-rbind(O,o)
        E<-E[-1,]
        E<-matrix(E,nrow=length(E)/2)
        if(length(E)==0){
          return(O)
          break
        }
      }
      else{
        TOUSE<-matrix(TOUSE,nrow=length(TOUSE)/2)
        o<-TOUSE[1,]
        TOUSE<-TOUSE[-1,]
      }
      if(length(E)!=0){
        E<-matrix(E,nrow=length(E)/2)
        l1<-which(E[,1]==o[2])
        l2<-which(E[,2]==o[2])
        for(i in l2){
          E[i,]<-rev(E[i,])
        }
        l<-append(l1,l2)
        E_o<-matrix(E[l,],nrow=length(E[l,])/2)
        if(length(E_o)==0){
        }
        else{
          E<-E[-l,]
          for(i in c(1:(length(E_o)/2))){
            if(abs(C[o[1],E_o[i,2]])<thres){
              O<-rbind(O,rev(E_o[i,]))
            }
            else{
              O<-rbind(O,E_o[i,])
              TOUSE<-rbind(TOUSE,E_o[i,])
            }
          }
        }
      }
    }
    else{
      return(O)
      break
    }
  }
}



##Auxiliary Functions

#Takes as input a graph g, and vector of edge-coefficients lambda (for polytrees)
#Returns the weight matrix L
graphtolambda<-function(g,lambda){
  E<-get.edgelist(g)
  n<-dim(E)[1]
  Lambda<-matrix(0,n+1,n+1)
  I<-diag(n+1)
  for(i in c(1:n)){
    Lambda[E[i,2],E[i,1]]<-lambda[i]
  }
  L<-solve(I-Lambda)
  return(L)
}

#Takes as input a graph g, and vector of edge-coefficients lambda, the size of the graph k and the list of oriented edges E (for DAGS)
#Returns the weight matrix L
DAGgraphtolambda<-function(O,k,g,lambda){
  e<-dim(E)[1]
  Lambda<-matrix(0,k,k)
  I<-diag(k)
  for(i in c(1:e)){
    Lambda[O[i,2],O[i,1]]<-lambda[i]
  }
  L<-solve(I-Lambda)
  return(L)
}


#Takes as input a size k.
#Returns both the adjacency matrix of a random polytree of size k and the adjacency matrix of its skeleton.
pruferwithskeleton <- function(k){
  if(k>2){
    P_orig <- ceiling(runif(k-2, min = 0, max = k))
    P <- P_orig
    V_orig <- 1:k
    V <- V_orig
    adj_matrix <- Matrix(matrix(numeric(k * k), ncol = k),sparse=TRUE)
    adj_matrixs <-Matrix(matrix(numeric(k * k), ncol = k),sparse=TRUE) 
    
    for(i in 1:(k-2)){
      complement <- setdiff(V, P)
      v_0 <- min(complement)
      V <- setdiff(V, v_0)
      adj_matrixs[which(V_orig == v_0), which(V_orig == P_orig[i])] <- adj_matrixs[which(V_orig == P_orig[i]), which(V_orig == v_0)]<-1
      
      m<-rbinom(1,1,0.5)
      if(m==0){
        adj_matrix[which(V_orig == v_0), which(V_orig == P_orig[i])] <- 1
      }
      else{
        adj_matrix[which(V_orig == P_orig[i]), which(V_orig == v_0)] <- 1
      }
      P <- P[2:length(P)]
    }
    adj_matrixs[which(V_orig == V[1]), which(V_orig == V[2])]<-adj_matrixs[which(V_orig == V[2]), which(V_orig == V[1])]<-1
    m<-rbinom(1,1,0.5)
    if(m==0){
      adj_matrix[which(V_orig == V[1]), which(V_orig == V[2])] <- 1
    }
    else{
      adj_matrix[which(V_orig == V[2]), which(V_orig == V[1])] <- 1
    }
  }
  else{
    adj_matrixs<-matrix(c(0,1,1,0),nrow=2,byrow=TRUE)
    m<-rbinom(1,1,0.5)
    if(m==0){
      adj_matrix<-matrix(c(0,1,0,0),nrow=2,byrow=TRUE)
    }
    else{
      adj_matrix<-matrix(c(0,0,1,0),nrow=2,byrow=TRUE)
    }
  }
  ret<-list(Directed=adj_matrix,Skeleton=adj_matrixs)
  return(ret)
}

#Takes as input a list of oriented edges.
#Returns an ajacency matrix.
edgelist_toadjmatrix<-function(O){
  n<-dim(O)[1]
  #I<-matrix(0,nrow=n+1,ncol=n+1)
  I<-Matrix(nrow=(n+1),ncol=(n+1),data=0,sparse=TRUE)
  for(i in c(1:n)){
    I[O[i,1],O[i,2]]<-1
  }
  return(I)
}

#Takes as input a data matrix x, and the position of three random variables.
#Returns E[x_ix_jx_k]
tmom<-function(x,i,j,k){
  n<-dim(x)[1]
  c<-c(i,j,k)
  x<-x[,c]
  M<-matrix(0,nrow=3,ncol=1)
  for(i in c(1:3)){
    M[i]<-mean(x[,i])
  }
  t<-0
  for (i in c(1:n)){
    t<-t+(1/n)*(x[i,1]-M[1])*(x[i,2]-M[2])*(x[i,3]-M[3])
  }
  return(t)
}

#Takes as input a data matrix x, and the position of four random variables.
#Returns E[x_ix_jx_kx_l]
fmom<-function(x,i,j,k,l){
  n<-dim(x)[1]
  c<-c(i,j,k,l)
  x<-x[,c]
  M<-matrix(0,nrow=4,ncol=1)
  for(i in c(1:4)){
    M[i]<-mean(x[,i])
  }
  t<-0
  for (i in c(1:n)){
    t<-t+(1/n)*(x[i,1]-M[1])*(x[i,2]-M[2])*(x[i,3]-M[3])*(x[i,4]-M[4])
  }
  return(t)
}

#Takes as input the adjacency matrix of the true (It) and the estimated graph (Ie).
#Returns the number of correctly/wrongly oriented, missing/extra edges and the false discovery rate between the true and estimated graphs.
diff_stat<-function(Ie,It){
  Ae<-Ie+t(Ie)
  At<-It+t(It)
  A<-Ae-At
  extra<-sum(A>0)/2
  missing<-sum(A<0)/2
  I<-Ie-It
  wrong<-sum(I>0)-extra
  correct<-sum(Ie>0)-extra-wrong
  if((wrong+extra+correct)==0){
    fdr=0
  }
  else{
    fdr<-(wrong+extra)/(wrong+extra+correct)
  }
  L<-list(CORRECT=correct,WRONG=wrong,MISSING=missing,EXTRA=extra,FDR=fdr)
  return(L)
}
