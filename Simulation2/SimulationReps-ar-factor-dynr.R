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

library(dynr)

for(rep in 1:200){
  set.seed(1004+rep*100)
  
  nT <- 100
  x0 <- rnorm(1,0,1)
  x <- c()
  x[1] <- x0
  mu1=0
  mu2=0
  
  ar1=0.1
  ar2=0.8
  
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
  
  data1 <- dynr.data(simdata, id="id", time="Time", observed=c("y.1","y.2","y.3"))
  
  meas <- prep.measurement(
    values.load=matrix(c(1,0,0,0,0,0),ncol=2), 
    params.load=matrix(c("fixed","l1","l2","fixed","fixed","fixed"),ncol=2,nrow=3),
    state.names=c("x","beta"),
    obs.names=c("y.1","y.2","y.3")) 
  
  initial <- prep.initial(
    values.inistate=c(1,.5), #a0
    params.inistate=c("mu", "beta0"),
    values.inicov=diag(c(1,0)), #P0
    params.inicov=diag(c("v1","fixed"))
  )
  
  mdcov <- prep.noise(
    values.latent=diag(c(1,0.1)), 
    params.latent=matrix(c("v1","fixed",
                           "fixed","v2"),ncol=2,byrow=T),
    values.observed=diag(c(1,1,1),3), 
    params.observed=diag(c("m1","m1","m1"),3))
  
  
  sample_formula = list(
    x ~ mu+beta*(x-mu),
    beta ~ beta
  )
  sample_dynamics = prep.formulaDynamics(formula = sample_formula,
                                         startval = c(mu = 1))
  
  VAR1model1 <- dynr.model(dynamics=sample_dynamics, measurement=meas, noise=mdcov,initial=initial, data=data1)
  VAR1model1@ub['beta0'] <- 1
  VAR1model1@lb['beta0'] <- -1
  VAR1model1@ub['mu'] <- max(y)
  VAR1model1@lb['mu'] <- min(y)
  
  res1 <- try(dynr.cook(VAR1model1,verbose=T,debug=T))
  if(class(res1)!="try-error"){
    summary(res1)
    round(coef(res1),3)
    
    # plot(x,type="b")
    # lines(res1@eta_smooth_final[1,],col=2,lty=3,lwd=2,type="l")
    # plot(res1@eta_smooth_final[2,],col=2,type="l",
    #      xlab="Observations",ylab=bquote("Smoothed"~ mu~" trajectory"))
    # lines(c(rep(ar1,nT/2),rep(ar2,nT/2)),lty=3,lwd=2)

    taste_VAR1<- dynr.taste(dynrModel= VAR1model1, dynrCook=res1, conf.level=.95, debug_flag = TRUE)
    shocks = data.frame(ID=taste_VAR1$id, Time=taste_VAR1$time, 
                        chi.jnt = as.numeric(taste_VAR1$chi.jnt.shock),
                        chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                        chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                        t.inn.x =as.numeric(taste_VAR1$t.inn.shock[1,]),#whether t-values are sig for V1
                        t.inn.ar =as.numeric(taste_VAR1$t.inn.shock[2,]) #whether t-values are sig for V2
    ) 
    
   save.image(paste0("Res/Res",rep,".RData"))
  }
 } 
  
