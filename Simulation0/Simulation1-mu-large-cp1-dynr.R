rm(list=ls())
library(dynr)

# 200 MC reps

for(rep in 1:200){
  
  
  cat(paste0("Starting replication ", rep,"\n"))
  
  set.seed(1004+rep*100)
  
  #### Data generation ####
  nT <- 100
  x0 <- rnorm(1,0,1)
  x <- c()
  x[1] <- x0
  mu1=0
  mu2=2
  
  ar1=0.3
  ar2=0.3
  
  for(i in 2:(nT/2)){
    x[i]<-mu1+ar1*(x[i-1]-mu1)+rnorm(1,0,1)
  }
  for(i in (nT/2+1):100){
    x[i]<-mu2+ar2*(x[i-1]-mu2)+rnorm(1,0,1)
  }
  
  simdata <- data.frame(x=x,id=1,Time=1:nT)
  save(simdata,file=paste0("SimRes/Sim1MuLargeCP1/AR03/Data/data_sim1_mu_large_rep",rep,".RData"))
  
  
  
  ###### dynr - SSM model fitting ###### 
  
  data1 <- dynr.data(simdata, id="id", time="Time", observed=c("x"))
  
  meas <- prep.measurement(
    values.load=matrix(c(1,0),ncol=2), 
    params.load=matrix("fixed",ncol=2,nrow=1),
    state.names=c("x","mu"),
    obs.names=c("x")) 
  
  initial <- prep.initial(
    values.inistate=c(0,0.8), #a0
    params.inistate=c("fixed", "mu0"),
    values.inicov=diag(c(1,0)), #P0
    params.inicov=diag(c("v1","fixed"))
  )
  
  mdcov <- prep.noise(
    values.latent=diag(c(1,0.1)), 
    params.latent=matrix(c("v1","fixed",
                           "fixed","v2"),ncol=2,byrow=T),
    values.observed=1e-4, 
    params.observed=matrix("fixed",ncol=1))
  
  
  sample_formula = list(
    x ~ mu+beta*(x-mu),
    mu ~ mu
  )
  sample_dynamics = prep.formulaDynamics(formula = sample_formula,
                                         startval = c(beta = .8))
  
  VAR1model1 <- dynr.model(dynamics=sample_dynamics, measurement=meas, noise=mdcov,initial=initial, data=data1)
  res1 <- dynr.cook(VAR1model1,verbose=F,debug=T)
  summary(res1)
  round(coef(res1),3)
  
  
  ### For plotting smoothed estimates for the TVP from FIS ###
  # plot(res1@eta_smooth_final[2,],col=2,type="l",
  #      xlab="Observations",ylab=bquote("Smoothed"~ mu~" trajectory"))
  # lines(c(rep(mu1,nT/2),rep(mu2,nT/2)),lty=3,lwd=2)
  
  
  ###### Deriving test statistics for outliers ###### 
  
  taste_VAR1<- dynr.taste(dynrModel= VAR1model1, dynrCook=res1, conf.level=.95, debug_flag = TRUE)
  
  
  shocks = data.frame(ID=taste_VAR1$id, Time=taste_VAR1$time, 
                      chi.jnt = as.numeric(taste_VAR1$chi.jnt.shock),
                      chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                      chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                      t.inn.x =as.numeric(taste_VAR1$t.inn.shock[1,]),
                      t.inn.mu =as.numeric(taste_VAR1$t.inn.shock[2,]) 
  ) 
  
  
  ### For plotting the outlier diagnostics ###
  
  # taste_plot <- autoplot(taste_VAR1)
  # library(gridExtra)
  # grid.arrange(taste_plot$`1`$chi_jnt,
  #              taste_plot$`1`$chi_inn,
  #              taste_plot$`1`$chi_add,
  #              nrow=3)
  # grid.arrange(taste_plot$`1`$`t_[x]`,
  #              taste_plot$`1`$`t_[mu]`,
  #              nrow=2)
  
  save(res1,taste_VAR1,shocks,file=paste0("SimRes/Sim1MuLargeCP1/AR03/dynr_sim1_mu_large_rep",rep,".RData"))

}

q("no")