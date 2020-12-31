inter.adj<-function(main.adj){
  n1=dim(main.adj)[1]
  adj2=main.adj
  diag(adj2)=1
  p.temp=(n1-1)*n1/2
  loc=location(n1)
  adj_inter=matrix(0,p.temp,p.temp)
  temp1=temp2=numeric()

  for(i in 1:(p.temp)){
    for(j in 1:p.temp){
      temp1=which(loc==i,arr.ind = TRUE)
      temp2=which(loc==j,arr.ind = TRUE)
      intersec=intersect(temp1,temp2)
      ii=setdiff(temp1,intersec)
      jj=setdiff(temp2,intersec)
      if(length(intersec)>0){
        tt=length(intersect(adj2[ii,],adj2[jj]))
        adj_inter[i,j]=ifelse(tt>0,1,0)
      }

    }
  }
  return(adj_inter)
}