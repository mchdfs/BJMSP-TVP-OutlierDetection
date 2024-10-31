rm(list=ls())

args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
  rep=0
  ##supply default values
}else{
  for(i1 in 1:length(args)){
    eval(parse(text=args[[i1]]))
  }
}

source("GAMfunctions_AR.R")

nT <- 100


for(rep in 2:200){
  
  cat(paste0("Starting replication ", rep,"\n"))
  
  load(file=paste0("SimRes/Sim3ARLargeCP1/0007/Data/data_sim3_AR_large_rep",rep,".RData"))
  
  data1 <- data.frame(x=c(simdata$x[-1],NA),lagx=simdata$x,Time=1:nT)
  
  temp5_once_AIC <- FindCP_gam(data1,threshold=-5)
  temp5_once_BIC <- FindCP_gam(data1,threshold=-5,IC="BIC")
  
  temp5AIC <- DynamicUpdate(data1)
  temp5BIC <- DynamicUpdate(data1,IC="BIC")
  
  save(temp5_once_AIC,temp5_once_BIC,temp5AIC,temp5BIC,file=paste0("SimRes/Sim3ARLargeCP1/new_mgcv_sim3_AR_large_rep",rep,".RData"))
  
}

q("no")