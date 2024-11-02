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

library(tidyr)
library(dplyr)

N=10
nT <- 100

for(rep in 301:350){
  
  set.seed(1004+rep*100)
  simdata <- c()
  for(n in 1:N){
    x0 <- rnorm(1,0,1)
    x <- c()
    x[1] <- x0
    mu1=0
    mu2=0
    
    ar1=0.1
    ar2=0.8
    breakpoint <- sample(15:(nT-15),1)
    
    for(i in 2:breakpoint){
      x[i]<-mu1+ar1*(x[i-1]-mu1)+rnorm(1,0,1)
    }
    for(i in (breakpoint+1):nT){
      x[i]<-mu2+ar2*(x[i-1]-mu2)+rnorm(1,0,1)
      
    }
    simdatai <- data.frame(x=x,id=n,Time=1:nT,breakpoint=breakpoint)
    simdata <- rbind(simdata,simdatai)
  }

  simdata$y <- simdata$x + rnorm(length(simdata$x),0,1)

###### dynr ------

library(dynr)
data1 <- dynr.data(simdata, id="id", time="Time", observed=c("y"))

meas <- prep.measurement(
  values.load=matrix(c(1,0),ncol=2), 
  params.load=matrix("fixed",ncol=2,nrow=1),
  state.names=c("x","beta"),
  obs.names=c("y")) 

initial <- prep.initial(
  values.inistate=c(0,0.8), #a0
  params.inistate=c("fixed", "beta0"),
  values.inicov=diag(c(1,0)), #P0
  params.inicov=diag(c("v1","fixed"))
)

mdcov <- prep.noise(
  values.latent=diag(c(1,0.1)), 
  params.latent=matrix(c("v1","fixed",
                         "fixed","v2"),ncol=2,byrow=T),
  values.observed=1e-4, 
  params.observed=matrix("vmeas",ncol=1))



sample_formula = list(
  x ~ mu+beta*(x-mu),
  beta ~ beta
)
sample_dynamics = prep.formulaDynamics(formula = sample_formula,
                                       startval = c(mu = 1))



VAR1model1 <- dynr.model(dynamics=sample_dynamics, measurement=meas, noise=mdcov,initial=initial, data=data1)
VAR1model1@ub['beta0'] <- 1
VAR1model1@lb['beta0'] <- -1
VAR1model1@ub['mu'] <- max(simdata$y)
VAR1model1@lb['mu'] <- min(simdata$y)

res1 <- try(dynr.cook(VAR1model1,verbose=T,debug=T))
  
  if(class(res1)!="try-error"){
    taste_VAR1<- try(dynr.taste(dynrModel= VAR1model1, dynrCook=res1, conf.level=.95, debug_flag = TRUE),silent=TRUE)
    if(class(taste_VAR1)!="try-error"){
      shocks = data.frame(ID=taste_VAR1$id, Time=taste_VAR1$time, 
                          chi.jnt = as.numeric(taste_VAR1$chi.jnt.shock),
                          chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                          chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                          t.inn.x =as.numeric(taste_VAR1$t.inn.shock[1,]),#whether t-values are sig for V1
                          t.inn.ar =as.numeric(taste_VAR1$t.inn.shock[2,]) #whether t-values are sig for V2
      ) 
      
      shocks %>% group_by(ID) %>% summarise(which(t.inn.ar>0)) -> TVPoutliers
      save.image(paste0("Res/Res",rep,".RData"))
    }
  }
}
