##################################################################################
####################Network Generation File#######################################
##################################################################################
library(plyr)
library(dplyr)
library(matrixStats)
library(R.matlab)
library(tidyr)
library(purrr)
numNetworks=300 #number of network types

ARRAYID=Sys.getenv('SLURM_ARRAY_TASK_ID') #from cluster batch job

oo=1
netParams<-read.csv("lhsNetParams.csv") #load network type parameters file

N<- 100000 #total population
HH_Size = 5 #household size

mat<-matrix(0, ncol=N/HH_Size, nrow=N/HH_Size)
e<-upper.tri(mat, diag = FALSE) 
pair_df<-as.data.frame(which(e==1, arr.ind = T)) #all possible HH connections
colnames(pair_df)<-c("a","b")
rm(e)
rm(mat)

#Gaussian Connectivity Kernel calculation adapted from Lang et al. 2018
connectivity_kernel <- function(distance, rho, n, sigma) {
  gKernel = exp((-distance^2)/(2*sigma^2))/(2*pi*sigma^2)
  Kernel= (n-4)*gKernel/rho #fixed network density (rho); individuals/area
  K= rowMins(as.matrix(data.frame(Kernel,1)))
  return(K)
}

d <- sqrt(N/HH_Size)

#these will be used to record degree, connection radius of actual genearted networks 
#as these will be slightly different from parameters used to generate networks
NetworkMetrics<-matrix(NA,nrow=dim(netParams)[1], ncol=5)
colnames(NetworkMetrics)<-c("ID","LHS_sigma", "LHS_deg","sigma", "deg")
NetworkMetrics[,1]<-1:dim(netParams)[1]
NetworkMetrics[,2]<-netParams[,2]
NetworkMetrics[,3]<-netParams[,3]

for (i in 1:dim(netParams)[1])#)
{
n<-netParams$deg[i]
sigma<-netParams$sigma[i]
    
## Add coordinates to pair dataframe
x <- runif(N/HH_Size, min = 0, max = d) # 1 by 10,000 max 100, min 0
y <- runif(N/HH_Size, min = 0, max = d)
    
df <- data.frame(id = 1:(N/HH_Size),
               x = x,
               y = y)
    
#save locations of all households
write.table(data.frame("0"=rownames(df),df),paste("Network", oo, "coordinates",ARRAYID, ".csv", sep=""), row.names=FALSE, sep=",")
  
## Extract the x and y coords for each pair
## for distance calculation
xa <- x[pair_df$a]
xb <- x[pair_df$b]
    
ya <- y[pair_df$a]
yb <- y[pair_df$b]

l=1
n1=rep(n[l],dim(pair_df)[1])
 
k=1 
sigma1 = sigma[k]

dd <- sqrt((xa-xb)^2 + (ya-yb)^2)
pair_df$dist <- dd
pair_df$na<-n1[pair_df$a] #average degree for house a
pair_df$nb<-n1[pair_df$b] #average degree for house b

area = sqrt(N/HH_Size)*sqrt(N/HH_Size)
rho= N/area #fixed

#calculate connectivity kernel for all HH pairs
pair_df$k<-connectivity_kernel(pair_df$dist, rho, n1, sigma1)

#total connects per HH pair based on K
draw_connects<-rbinom(length(pair_df$k),(HH_Size*HH_Size),pair_df$k)

#subset only connected houses
connectedHHs<-pair_df[draw_connects>0,]

#add total connections to connected houses, i.e., there can be a max of 25 connections between 2 houses
connectedHHs$totConnects<-draw_connects[draw_connects>0]
HHcombos<-matrix(c(1,1, 1,2, 1,3, 1,4, 1,5, 2,1, 3,1, 4,1, 5,1, 2,2, 2,3, 2,4, 2,
                   5, 5,2, 4,2, 3,2, 3,3, 3,4, 3,5, 5,3, 4,3, 4,4, 4,5, 5,4, 5,5), ncol=2, byrow=T)

network<-c()
for (nrConn in 1:12) #hardcoded half way 
{
  #print(nrConn)
  Conn<-connectedHHs[connectedHHs$totConnects==nrConn,]
  if (dim(Conn)[1]>0){
    nrSample<-matrix(sample(1:25,dim(Conn)[1]*nrConn, replace=T), ncol=nrConn) #sample two times
    
    #check and remove duplicates e.g. if the same connection was selected twice within same HH pair
    check<-apply(nrSample, 1, function(x) length(unique(x)))-nrConn
    while (sum(check<0)>0) #if negative then duplicate
    {
      nrSample[which(check<0),]<-matrix(sample(1:25,sum(check<0)*nrConn, replace=T), ncol=nrConn) #replace second one with new sample
      check<-apply(nrSample, 1, function(x) length(unique(x)))-nrConn #recheck
    }
    
    Conn<-cbind(Conn[,1:3],nrSample)
    
    #get PIDs from first and second interhouse connections 
    PIDConn1<-HHcombos[unlist(Conn[,seq(4,nrConn+3,1)]),]#find corresponding PIDs from combo lookup table
    a<-(Conn[,1]-1)*5+PIDConn1[,1] #convert to PID 
    b<-(Conn[,2]-1)*5+PIDConn1[,2] #convert to PID
    dist<-rep(Conn[,3], nrConn) #don't need each, HHs alternate
    network<-rbind(network,cbind(a,b, dist))
    
  }
}

network1=c()
for (nrConn in 13:24) #hardcoded other half
{
  #print(nrConn)
  Conn<-connectedHHs[connectedHHs$totConnects==nrConn,]
  if (dim(Conn)[1]>0){
    nrSample<-matrix(sample(1:25,dim(Conn)[1]*(25-nrConn), replace=T), ncol=(25-nrConn)) #select non-sampled
    #check and remove duplicates e.g. if the same connection was selected twice within same HH pair
    check<-apply(nrSample, 1, function(x) length(unique(x)))-(25-nrConn)
    while (sum(check<0)>0) #if negative then duplicate
    {
      nrSample[which(check<0),]<-matrix(sample(1:25,sum(check<0)*(25-nrConn), replace=T), ncol=(25-nrConn)) #replace second one with new sample
      check<-apply(nrSample, 1, function(x) length(unique(x)))-(25-nrConn) #recheck
    }
    
    fullSet<- matrix(rep(1:25 , dim(Conn)[1]), nrow=dim(Conn)[1], byrow=T)
    rowindx<-matrix(rep(1:dim(Conn)[1],each=(25-nrConn)), nrow=dim(Conn)[1], byrow=T)
    all<-c()
    for (j in 1:dim(rowindx)[1])
    {
      selectedConns<-c()
      fullSet[rowindx[j,],nrSample[j,]]<-0
      
    }
    p<-c(t(fullSet)) #restructure
    g<-p[p!=0] #remove zeros
    selectedConns<-matrix(g, ncol=nrConn, byrow=T)
    Conn<-cbind(Conn[,1:3], selectedConns)
    
    #get PIDs from first and second interhouse connections 
    PIDConn1<-HHcombos[unlist(Conn[,seq(4,nrConn+3,1)]),]#find corresponding PIDs from combo lookup table
    a<-(Conn[,1]-1)*5+PIDConn1[,1] #convert to PID 
    b<-(Conn[,2]-1)*5+PIDConn1[,2] #convert to PID
    dist<-rep(Conn[,3], nrConn) #don't need each, HHs alternate
    network1<-rbind(network1,cbind(a,b, dist))
  }
}


#for 25 connections
nrConn=25
Conn<-connectedHHs[connectedHHs$totConnects==25,]
network2=c()
if (dim(Conn)[1]>0){
Conn<-cbind(Conn[,1:3], matrix(rep(1:25 , dim(Conn)[1]), nrow=dim(Conn)[1], byrow=T))
PIDConn1<-HHcombos[unlist(Conn[,seq(4,nrConn+3,1)]),]#find corresponding PIDs from combo lookup table
a<-(Conn[,1]-1)*5+PIDConn1[,1] #convert to PID 
b<-(Conn[,2]-1)*5+PIDConn1[,2] #convert to PID
dist<-rep(Conn[,3], nrConn) #don't need each, HHs alternate
network2<-cbind(a,b, dist)
}


#add HH
HH_ID<-N/HH_Size
HH_ID<- 1:HH_ID

network_df<-rbind(network,network1,network2)  
nodes <- data.frame(x =c(network_df[,1],network_df[,2]))
degree <- nodes %>% group_by(x) %>%  summarise(degree = n()) 

colnames(network_df)<-c("conn1","conn2","dist")
network_df<-data.frame(network_df)

network_df2<-network_df[,c(2,1,3)]

colnames(network_df2)<-c("conn1","conn2","dist")

NetMat<-data.frame(rbind(network_df,network_df2 ))

NetMat<-NetMat[order(NetMat$conn1),]

#######################################################################################
#Save Network, ########################################################################
#######################################################################################



write.csv(NetMat,paste("Network", oo,"array",ARRAYID, ".csv", sep=""))


oo=oo+1


NetworkMetrics[i,4]<-mean(network_df[,3])
NetworkMetrics[i,5]<-mean(degree$degree)


}

#save avearge dergee and connection radius
write.csv(NetworkMetrics,paste("NetworkMetrics",ARRAYID, ".csv", sep=""))

