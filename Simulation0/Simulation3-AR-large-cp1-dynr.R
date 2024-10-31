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
for(rep in 21:200){
set.seed(1004+rep*100)

nT <- 100
x0 <- rnorm(1,0,1)
x <- c()
x[1] <- x0
mu1=0
mu2=0

ar1=0.8
ar2=0.1

for(i in 2:(nT/2)){
  x[i]<-mu1+ar1*(x[i-1]-mu1)+rnorm(1,0,1)
}
for(i in (nT/2+1):100){
  x[i]<-mu2+ar2*(x[i-1]-mu2)+rnorm(1,0,1)
}

simdata <- data.frame(x=x,id=1,Time=1:nT)
save(simdata,file=paste0("SimRes/Sim3ARLargeCP1/Data/data_sim3_AR_large_rep",rep,".RData"))


# plot(simdata$x,type="l",xlab="Time")
# abline(v=51,lty=2)



###### dynr ------

library(dynr)
data1 <- dynr.data(simdata, id="id", time="Time", observed=c("x"))

meas <- prep.measurement(
  values.load=matrix(c(1,0),ncol=2), 
  params.load=matrix("fixed",ncol=2,nrow=1),
  state.names=c("x","beta"),
  obs.names=c("x")) 

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
  params.observed=matrix("fixed",ncol=1))


sample_formula = list(
  x ~ mu+beta*(x-mu),
  beta ~ beta
)
sample_dynamics = prep.formulaDynamics(formula = sample_formula,
                                       startval = c(mu = .8))

VAR1model1 <- dynr.model(dynamics=sample_dynamics, measurement=meas, noise=mdcov,initial=initial, data=data1,outfile = paste0("cfile",rep,".c"))
res1 <- dynr.cook(VAR1model1,verbose=T,debug=T)
summary(res1)
round(coef(res1),3)


# plot(x,type="b")
# lines(res1@eta_smooth_final[1,],col=2,lty=3,lwd=2,type="l")
# plot(res1@eta_smooth_final[2,],col=2,type="l",
#      xlab="Observations",ylab=bquote("Smoothed"~ beta~" trajectory"))
# lines(c(rep(ar1,nT/2),rep(ar2,nT/2)),lty=3,lwd=2)


taste_VAR1<- dynr.taste(dynrModel= VAR1model1, dynrCook=res1, conf.level=.95, debug_flag = TRUE)
shocks = data.frame(ID=taste_VAR1$id, Time=taste_VAR1$time, 
                    chi.jnt = as.numeric(taste_VAR1$chi.jnt.shock),
                    chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                    chi.inn = as.numeric(taste_VAR1$chi.inn.shock),
                    t.inn.x =as.numeric(taste_VAR1$t.inn.shock[1,]),#whether t-values are sig for V1
                    t.inn.ar =as.numeric(taste_VAR1$t.inn.shock[2,]) #whether t-values are sig for V2
) 



# taste_plot <- autoplot(taste_VAR1)
# library(gridExtra)
# grid.arrange(taste_plot$`1`$chi_jnt,
#              taste_plot$`1`$chi_inn,
#              taste_plot$`1`$chi_add,
#              nrow=3)
# grid.arrange(taste_plot$`1`$`t_[x]`,
#              taste_plot$`1`$`t_[beta]`,
#              nrow=2)

save(res1,taste_VAR1,shocks,file=paste0("SimRes/Sim3ARLargeCP1/dynr_sim3_AR_large_rep",rep,".RData"))

}

q("no")
