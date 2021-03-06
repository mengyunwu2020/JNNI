% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GANAIr.path.R
\name{GANAIr.path}
\alias{GANAIr.path}
\title{fit a Bayesian model for the integration of nodes, interaction terms and networks}
\usage{
GANAIr.path(X, pathway, adjlist, t, s_init, a_init)
}
\arguments{
\item{X}{a matrix of predictor variable}

\item{pathway}{a vector of pathway information, for instance, value 1 represent the corresponding node location belongs to pathway 1}

\item{adjlist}{a list of Laplace information, each of which is a matrix capturing the relationship between node in the corresponding pathway}

\item{t}{a vector of the response variable}

\item{s_init}{a vector of tunning parameter related to the variance}

\item{a_init}{a vector of tunning parameter related to the initial selecting probability}
}
\value{
a list, each element corresponding to its tunning parameters contains values  as follows:
\item{mean_w_red}{The coefficients vector}
\item{beta_interaction}{The coefficients of interaction terms}
\item{beta_main}{The coefficients of main effects}
\item{BIC}{Bayesian Information Criterion}
\item{df}{The number of nonzero coefficients}
\item{index_path}{Selected networks' number}
}
\description{
One of the main functions in the JNNI package. Fits a path of JNNI models over different values of the tunning parameters.
}
\examples{
set.seed(11)
u=runif(4,0.8,1.2)
u=round(u,2)
rate=round(1/sqrt(5),2)
n=300
p=1000
rho=0.6
a1=c(u[1],rate*rep(u[1],5))
a2=c(u[2],rate*rep(u[2],5))
a3=c(u[3],rate*rep(u[3],5))
b1=rate*rep(u[1],5)
b2=rate*rep(u[2],5)
b3=rate*rep(u[3],2)
c1=rate*rep(u[4],5)
M=100
set.seed(1)
X=matrix(NA,n,p)
p0=p/M
for(i in 1:n){
trans=rnorm(M)
for(j in 1:M){
start = (j-1)*p0+1
X[i,start]= trans[j]
xx=matrix(rnorm((p0-1), mean =(trans[j]*rho), sd=sqrt(1-rho^2)),ncol=1)
X[i,(start+1):(start+(p0-1))]=xx
}
}
x=list()
for(i in 1:M){
x[[i]]=X[,(1+(i-1)*p0):(p0*i)]
}

Mm=M+M*(M-1)/2
Phix=list()
for(m in 1:M){
r=0
p_hat=p0+p0*(p0-1)/2
temp1=matrix(NA,n,p_hat)
temp1[,1:p0]=x[[m]]
for(i in 1:(p0-1)){
for(j in (1+i):p0){
r=r+1
temp1[,p0+r]=(x[[m]][,i])*(x[[m]][,j])
}
}
Phix[[m]]=temp1
}
cp=0
R1=matrix(0,M,M)
for(k in 1:(M-1)){
for(m in (k+1):M){
cp=cp+1
R1[k,m]=cp
temp1=matrix(NA,n,(p0*p0))
r=1
for(i in 1:p0){
temp1[,r:(r+p0-1)]=(x[[m]][,1:p0])*(x[[k]][,i])
r=r+p0
}
Phix[[M+cp]]=temp1
}
}

real_w=list()
for(i in 1:M)
real_w[[i]]=rep(0,times=p0+(p0-1)*p0/2)
for(i in (1+M):Mm){
indexp1=which(R1==(i-M),arr.ind = TRUE)[1]
indexp2=which(R1==(i-M),arr.ind = TRUE)[2]
real_w[[i]]=rep(0,times=p0*p0)
}
real_w[[1]][1:length(a1)]=a1
real_w[[1]][(p0+1):(p0+length(b1))]=b1
real_w[[2]][1:length(a2)]=a2
real_w[[2]][(p0+1):(p0+length(b2))]=b2
real_w[[3]][1:length(a3)]=a3
real_w[[3]][(p0+1):(p0+length(b3))]=b3
real_w[[M+1]][1:length(c1)]=c1


real_y=0
for(m in c(1:3,(M+1))){
real_y=Phix[[m]]\%*\%real_w[[m]]+real_y
}
t = real_y + rnorm(n)
adj=matrix(0,p0,p0)
for(i in 2:p0)
adj[1,i]=1
adj = adj +t(adj)
C=NULL
for(k in 1:M)
C[[k]] = adj


pathway=matrix(0,p,1)
p_each=p0
for (i in 1:M){
pathway[(p_each*(i-1)+1):(p_each*i),1]=i
}
s_init=2e-4
a_init=c(0.95,0.94)
temp_1<-GANAIr.path(X,pathway,C,t,s_init,a_init)
}
