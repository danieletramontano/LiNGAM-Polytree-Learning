source("Functions.R")

#Set parameter
p<-10 #number of nodes
n<-100 #sample size
d<-10 #max in degree (to be used only for general DAGs)
thres<-0.05 #Threshold for algorithms PT and TP.
k<-100 #number of repetitions 


#Matrix containing the (mean) statistics for the algorithms, in order Corrected oriented edges,
#wrongly oriented edges, missing edges, extra edges, False discovery rate and runtime (in seconds)
A<-matrix(0,nrow=5,ncol=6)
colnames(A)<-c("correct","wrong","extra","missing","FDR","time")
row.names(A)<-c("P","PT","TP","ANM","LiNGAM")


for(i in c(1:k)){
  #Polytree Generation, comment if you want to test the algorithm on a general DAG
  I<-pruferwithskeleton(p)
  Itestdirected<-I$Directed                                           #Adjacency matrix
  g<-graph_from_adjacency_matrix(Itestdirected)                       #graph
  
  
  #DAG generation, uncomment if you want to test the algorithm on a general DAG
  #g<-graph_from_graphnel(randDAG(p,d))                               #graph
  #Itestdirected<-as_adjacency_matrix(g)                              #Adjacency matrix
  
  E<-get.edgelist(g)              
  e<-dim(E)[1]                                                        #Number of edges
  E<-matrix(as(E,"numeric"),ncol=2)                                   #Edge list
  
  
  epsilon<-matrix(numeric(),nrow=n,ncol=p)                            #Error vector
  for(l in c(1:p)){
    a<-runif(1,-10,-1)
    s<-runif(1,1,10)
    epsilon[,l]<-runif(n,a,s)-0.5*(a+s)
  }
  
  lambda<-numeric()                                                   #Vector of coefficients
  for(l in c(1:e)){
    m<-rbinom(1,1,0.5)
    if(m==0){
      lambda[l]<-runif(1,-1,-0.3)
    }
    else{
      lambda[l]<-runif(1,0.3,1)
    }
  }
  
  W<-DAGgraphtolambda(E,p,g,lambda)                                   #Coefficient matrix
  x<-t(W%*%t(epsilon))                                                #Data matrix
  
  #CHOW LIU
  start_time<-Sys.time()
  S<-cov(x)                                                           #Covariance matrix
  C<-cov2cor(S)                                                       #Correlation matrix
  
  w<--abs(as.vector(C[t(upper.tri(C))]))                              #Vector of weights 
  g_full<-make_full_graph(p)       
  m_new<-mst(g_full,w)                                                #skeleton 
  E_n<-get.edgelist(m_new)                                            #List of estimated unoriented edges
  end_time<-Sys.time()
  time<-difftime(end_time,start_time,unit="secs")                                          #Chow-Liu runtime
  
  #PAIRWISE
  start_time<-Sys.time()
  Ip<-Oneedgeordering_third_fourth_cum(x,S,E_n)                       #Pairwise Adjacency matrix 
  end_time<-Sys.time()
  A["P",-6]<-A["P",-6]+as(diff_stat(Ip,Itestdirected),"numeric")/k                #Pairwise Statistic vector
  A["P",6]<-A["P",6]+time+difftime(end_time,start_time,unit="secs")/k                                   #Pairwise Runtime
  
  #PT
  start_time<-Sys.time()
  V_n<-c(1:p)
  LIST<-triplets(E_n,C,V_n,thres)
  Ulist<-LIST$Ulist
  Olist<-LIST$Olist
  
  Olist1<-RecursiveOneEdge_third_fourth_cum(x,S,Ulist,Olist)
  Itp<-edgelist_toadjmatrix(Olist1)                                   #PT Adjacency matrix
  end_time<-Sys.time()
  A["PT",-6]<-A["PT",-6]+as(diff_stat(Itp,Itestdirected),"numeric")/k             #PT Statistic vector
  A["PT",6]<-A["PT",6]+time+difftime(end_time,start_time,unit="secs")/k                                  #PT runtime
  
  #TP
  start_time<-Sys.time()
  Olistpt<-tmomcor_third_fourth_cum(x,C,S,E_n,thres)
  Ipt<-edgelist_toadjmatrix(Olistpt)                                  #TP Adjacency matrix
  end_time<-Sys.time() 
  A["TP",-6]<-A["TP",-6]+as(diff_stat(Ipt,Itestdirected),"numeric")/k             #TP Statistic vector
  A["TP",6]<-A["TP",6]+time+difftime(end_time,start_time,unit="secs")/k                                  #TP runtime
  
  #ANM
  start_time<-Sys.time()
  Ianm<-Oneedgeordering_anm(x,S,E_n)                                    #ANM Adjacency matrix 
  end_time<-Sys.time()
  A["ANM",-6]<-A["ANM",-6]+as(diff_stat(Ianm,Itestdirected),"numeric")/k            #ANM Statistic vector
  A["ANM",6]<-A["ANM",6]+difftime(end_time,start_time,unit="secs")/k     
  
  #LINGAM
  start_time<-Sys.time()
  Ilingam<-((t(lingam(x)$Bpruned))!=0)*1                                 #LiNGAM Adjacency matrix
  end_time<-Sys.time()
  A["LiNGAM",-6]<-A["LiNGAM",-6]+as(diff_stat(Ilingam,Itestdirected),"numeric")/k   #LiNGAM STATISTICS
  A["LiNGAM",6]<-A["LiNGAM",6]+difftime(end_time,start_time,unit="secs")/k
}


