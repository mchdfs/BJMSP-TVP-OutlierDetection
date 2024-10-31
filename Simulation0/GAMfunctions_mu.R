library(mgcv)

# Adapted from Albers & Bringmann Supplementary, 2020
# Current mgcv version 1.9-1

nullModel_gam<-function(Data,Var,VarL){
  tt=1:dim(Data)[1]
  N=dim(Data)[1]
  model0 <- gam(Var ~ s(tt) + VarL,data=Data)
  
  AIC0<- model0$aic
  BIC0<- BIC(model0)
  
  return(list(AIC0=AIC0,BIC0=BIC0,model0=model0))}


CPModel_gam<-function(Data,Var,VarL){
  D<-c()
  tt=1:dim(Data)[1]
  N=dim(Data)[1]
  N_1<-N-1
  AIC1 <- rep(NA,N)
  BIC1 <- rep(NA,N)
  
  modelall <- list()
  for(i in 1:N){
    D <- c(rep(0,i),rep(1,N-i)) 
    
    model1=gam(Var ~ s(tt) + D + VarL,data=Data)
    
    AIC1[i] <- model1$aic
    BIC1[i] <- BIC(model1)
    
    modelall[[i]] <- model1
  }
  return(list(AIC1=AIC1,BIC1=BIC1,model1=modelall))
}


FindCP_gam<- function(Data,threshold=-5,IC="AIC"){
  model0 <- nullModel_gam(Data=Data,Var=Data$x,VarL=Data$lagx)
  AIC0<- model0$AIC0
  BIC0<- model0$BIC0
  
  CPmodel1 <- CPModel_gam(Data=Data,Var=Data$x,VarL=Data$lagx)
  AIC1 <- CPmodel1$AIC1
  BIC1 <- CPmodel1$BIC1
  if(IC=="AIC"){
    return(list(CPdetectAIC = any(AIC1-AIC0 < threshold), 
                CPlocAIC=which.min(AIC1-AIC0),
                model0=model0,
                model1=CPmodel1,
                AICdiff=AIC1-AIC0))
  }else{if(IC=="BIC"){
    return(list(CPdetectBIC = any(BIC1-BIC0 < threshold), 
                CPlocBIC=which.min(BIC1-BIC0),
                model0=model0,
                model1=CPmodel1,
                BICdiff=BIC1-BIC0))
  }
  }
}


DynamicUpdate <- function(Data,threshold=-5,IC="AIC"){
  chunks_to_run <- matrix(c(1,100),nrow=1)
  depth <- 1
  if(IC=="AIC"){
    CPlocAIC_all <- c()
    repeat{
      res_temp <- try(FindCP_gam(Data[chunks_to_run[depth,1]:chunks_to_run[depth,2],],threshold=threshold))
      print(depth)
      print(CPlocAIC_all)
      print(chunks_to_run)
      if(class(res_temp)!="try-error"&&res_temp$CPdetect){
        CPlocAIC_all <- c(CPlocAIC_all,chunks_to_run[depth,1]+res_temp$CPlocAIC-1)
        if(abs(res_temp$CPlocAIC)>10){chunks_to_run <- rbind(chunks_to_run,c(chunks_to_run[depth,1],chunks_to_run[depth,1]+res_temp$CPlocAIC-2))}
        if((length(res_temp$AICdiff)-res_temp$CPlocAIC)>10){chunks_to_run <- rbind(chunks_to_run,c(chunks_to_run[depth,1]+res_temp$CPlocAIC,chunks_to_run[depth,2]))}
      }
      if(depth == nrow(chunks_to_run)){break}
      depth <- depth+1
    }
    return(list(CPlocAIC_all=CPlocAIC_all,chunks_to_run=chunks_to_run))
  }
  if(IC=="BIC"){
    CPlocBIC_all <- c()
    repeat{
      res_temp <- try(FindCP_gam(Data[chunks_to_run[depth,1]:chunks_to_run[depth,2],],threshold=threshold,IC="BIC"))
      print(depth)
      print(CPlocBIC_all)
      print(chunks_to_run)
      if(class(res_temp)!="try-error"&&res_temp$CPdetect){
        CPlocBIC_all <- c(CPlocBIC_all,chunks_to_run[depth,1]+res_temp$CPlocBIC-1)
        if(abs(res_temp$CPlocBIC)>10){chunks_to_run <- rbind(chunks_to_run,c(chunks_to_run[depth,1],chunks_to_run[depth,1]+res_temp$CPlocBIC-2))}
        if((length(res_temp$BICdiff)-res_temp$CPlocBIC)>10){chunks_to_run <- rbind(chunks_to_run,c(chunks_to_run[depth,1]+res_temp$CPlocBIC,chunks_to_run[depth,2]))}
      }
      if(depth == nrow(chunks_to_run)){break}
      depth <- depth+1
    }
    return(list(CPlocBIC_all=CPlocBIC_all,chunks_to_run=chunks_to_run))
  }  
  
}
