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


for(rep in 1201:1300){
  set.seed(1004+rep*100)
  
  nT <- 100
  x0 <- rnorm(1,0,1)
  x <- c()
  x[1] <- x0
  mu1=0
  mu2=2
  
  ar1=0.5
  ar2=0.5
  
  for(i in 2:(nT/2)){
    x[i]<-mu1+ar1*(x[i-1]-mu1)+rnorm(1,0,1)
  }
  for(i in (nT/2+1):nT){
    x[i]<-mu2+ar2*(x[i-1]-mu2)+rnorm(1,0,1)
  }
  
  measloading <- c(1,0.8,0.4)
  
  y <- t(measloading%*%t(x)) + rnorm(length(x)*3,0,1)
  
  simdata <- data.frame(x=x,y=y,id=1,Time=1:nT)
  
 
  
  ###### dynr ------
  
  library(dynr)
  data1 <- dynr.data(simdata, id="id", time="Time", observed=c("y.1","y.2","y.3"))
  
  meas <- prep.measurement(
    values.load=matrix(c(1,0,0,0,0,0),ncol=2), 
    params.load=matrix(c("fixed","l1","l2","fixed","fixed","fixed"),ncol=2,nrow=3),
    state.names=c("x","mu"),
    obs.names=c("y.1","y.2","y.3")) 
  
  initial <- prep.initial(
    values.inistate=c(0,1), #a0
    params.inistate=c("fixed", "mu0"),
    values.inicov=diag(c(1,0)), #P0
    params.inicov=diag(c("v1","fixed"))
  )
  
  mdcov <- prep.noise(
    values.latent=diag(c(1,0.1)), 
    params.latent=matrix(c("v1","fixed",
                           "fixed","v2"),ncol=2,byrow=T),
    values.observed=diag(c(0.1,0.1,0.1),3), 
    params.observed=diag(c("m1","m2","m3"),3))
  
  
  sample_formula = list(
    x ~ mu+beta*(x-mu),
    mu ~ mu
  )
  sample_dynamics = prep.formulaDynamics(formula = sample_formula,
                                         startval = c(beta = .8))
  
  VAR1model1 <- dynr.model(dynamics=sample_dynamics, measurement=meas, noise=mdcov,initial=initial, data=data1)
  VAR1model1@ub['beta'] <- 1
  VAR1model1@lb['beta'] <- -1
  VAR1model1@ub['mu0'] <- max(y)
  VAR1model1@lb['mu0'] <- min(y)
  
  res1 <- try(dynr.cook(VAR1model1,verbose=T,debug=T),silent=TRUE)
 
  
  if(class(res1)!="try-error"){
    taste_VAR1<- try(dynr.taste(dynrModel= VAR1model1, dynrCook=res1, conf.level=.95, debug_flag = TRUE),silent=TRUE)
    
    if(class(taste_VAR1)!="try-error"){
      shocks = data.frame(ID=taste_VAR1$id, Time=taste_VAR1$time, 
                        chi.jnt = as.numeric(taste_VAR1$chi.jnt.shock),
                        chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                        t.inn.x =as.numeric(taste_VAR1$t.inn.shock[1,]),#whether t-values are sig for V1
                        t.inn.ar =as.numeric(taste_VAR1$t.inn.shock[2,]) #whether t-values are sig for V2
    
    ) 
    
    

    save.image(paste0("Res/Res",rep,".RData"))
    }
  }
}

q("no")
