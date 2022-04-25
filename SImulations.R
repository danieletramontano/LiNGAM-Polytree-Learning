source("AllTheFunctions_sparse.R")


#Set parameter
p<-10 #number of nodes
n<-20 #sample size
d<-10 #max in degree (to be used only for general DAGs)
thres<-0.05 #Threshold for algorithms PT and TP.
k<-10 #number of repetitions 

#Vector containing the (mean) statistics for the algorithms, in order Corrected oriented edges, wrongly oriented edges, missing edges, extra edges and False discovery rate
S_p<-S_pt<-S_tp<-matrix(0,nrow=1,ncol=5)
#Runtime
T_p<-T_pt<-T_tp<-0

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
  E_n<-get.edgelist(m_new)                                            #List of unoriented edges
  end_time<-Sys.time()
  time<-end_time-start_time                                           #Chow-Liu runtime
  
  #PAIRWISE
  start_time<-Sys.time()
  Ip<-Oneedgeordering_third_fourth_cum(x,S,E_n)                       #Adjacency matrix 
  end_time<-Sys.time()
  S_p<-S_p+as(diff_stat(Ip,Itestdirected),"numeric")/k                #Pairwise Statistic vector
  T_p<-time+(end_time-start_time)/k                                   #Pairwise Runtime
  
  #PT
  start_time<-Sys.time()
  V_n<-c(1:p)
  LIST<-triplets(E_n,C,V_n,thres)
  Ulist<-LIST$Ulist
  Olist<-LIST$Olist
  
  Olist1<-RecursiveOneEdge_third_fourth_cum(x,S,Ulist,Olist)
  Itp<-edgelist_toadjmatrix(Olist1)                                   #Adjacency matrix
  end_time<-Sys.time()
  S_pt<-S_pt+as(diff_stat(Itp,Itestdirected),"numeric")/k             #PT Statistic vector
  T_pt<-time+(end_time-start_time)/k                                  #PT runtime
  
  #TP
  start_time<-Sys.time()
  Olistpt<-tmomcor_third_fourth_cum(x,C,S,E_n,thres)
  Ipt<-edgelist_toadjmatrix(Olistpt)                                  #Adjacency matrix
  end_time<-Sys.time() 
  S_tp<-S_tp+as(diff_stat(Ipt,Itestdirected),"numeric")/k             #TP Statistic vector
  T_tp<-time+(end_time-start_time)/k                                  #TP runtime
}
