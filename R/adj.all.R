adj.all<-function(adj){
  t1=Sys.time()
  temp=list()
  temp[[1]]=adj
  temp[[2]]=inter.adj(adj)
  adj.all=as.matrix(Matrix::bdiag(temp))
  # t2=Sys.time()
  # print(t2-t1)

  # t1=Sys.time()
  # p1=dim(adj)[1]
  # p2=p1*(p1+1)/2
  # adj.all=matrix(0,p2,p2)
  # adj.all[1:p1,1:p1]=adj
  # adj.all[(p1+1):p2,(p1+1):p2]=inter.adj(adj)
  # t2=Sys.time()
  # print(t2-t1)
  return(adj.all)
}