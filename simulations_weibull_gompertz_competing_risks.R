# Testing similarity of parametric competing risks models for identifying potentially similar pathways in healthcare #
# https://doi.org/10.48550/arXiv.2401.04490
# Transition intensities based on the Gompertz distribution (state 1 and 2) and the Weibull distribution (state 3)
# Administrative censoring (here after 90 days)
# Simulation Study based on application example (Scenario 3 and Scenario 4 in the paper)
# Code uses a parallelization as the simulation is quite computational intensive
# Author: Kathrin Moellenhoff

library(alabama)
library(dplyr)
library(flexsurv)
library(foreach)
library(doParallel)

######
# simulate data for two groups (based on the transition intensities of the application example)
# specify n1 and n2 in advance
# specify thresholds for the test on similarity
# fix the level alpha (here 0.05)
# fix the number of bootstrap repetitions (here B=250)
# specify the end of study (here 90 days)
#####

n1 <- 300 #sample sizes
n2 <- 300

alpha <- 0.05 #significance level
B <- 250 #bootstrap repetitions
Nsim <- 1000 #number of simulations

# function for data generation
sim.csh.gompertzandweibull <- function(n,
                                       shape.gomp,  # Gompertz shape (vector for each hazard)
                                       rate.gomp,  # Gompertz rate/scale >0 (vector for each hazard)
                                       shape.weib, # Weibull shape
                                       scale.weib, # Weibull scale
                                       max.int=5000, #for call of the function uniroot, the interval is adjusted by extend.Int (so it can not be chosen too small)
                                       cens.time # for administrative censoring, no censoring for cens.time=NULL 
)
{
  
  #########################
  # event.times
  
  #cumulative all-cause hazard + y
  temp <- function(x,y, shape.gomp, rate.gomp, shape.weib, scale.weib)
  {return(sum(Hgompertz(x,shape=shape.gomp, rate=rate.gomp))+sum(Hweibull(x,shape=shape.weib, scale=scale.weib))+y)} # Hgompertz and Hweibull from flexsurv package for cumulative hazards
  
  event.times <- NULL
  i <- 1
  n.failure <- 0
  
  #generate event times
  while(i <= n)
  {
    u <- runif(1)
    res <- uniroot(temp, c(0,max.int), tol=0.000000001, y=log(1-u), shape.gomp=shape.gomp, rate.gomp=rate.gomp,shape.weib=shape.weib, scale.weib=scale.weib,extendInt = "upX") #oder extendInt="yes
    event.times[i] <- res$root
    i <- i+1
  }
  
  ##################################################
  ##failure cause
  
  #total number of competing risks 
  n.comp.risks<-length(shape.gomp)+length(shape.weib)
  
  event.probs = matrix(NA, nrow=n, ncol = n.comp.risks)
  for(i in 1:length(shape.gomp)){
    event.probs[,i] = hgompertz(event.times,shape=shape.gomp[i],rate=rate.gomp[i])
  }
  for(i in 1:length(shape.weib)){
    event.probs[,i+length(shape.gomp)] = hweibull(event.times,shape=shape.weib[i],scale=scale.weib[i])
  }
  
  #binomial experiment to determine failure cause (f.cause)
  f.cause = matrix(NA, nrow=n.comp.risks, ncol = n)
  for(i in 1:n){
    f.cause[,i] <- rmultinom(1,
                             size=1,
                             prob = event.probs[i,] / rowSums(event.probs)[i]
    )
  }
  
  f.cause =  as.numeric(1:n.comp.risks %*% f.cause)
  
  # administrative censoring at time cens.time
  if(is.null(cens.time)){
    obs.times <- event.times
    obs.cause <- f.cause
  } 
  else{
    obs.times <- pmin(event.times, cens.time)
    obs.cause <- c(event.times <= cens.time)*f.cause
  }
  
  ##################################################
  
  crd <- data.frame(id=1:n,
                    from = rep(0,n),
                    to = obs.cause, 
                    time = obs.times 
  )
  
  return(crd)
}  

# this is needed for parallelizing 
numCores <- detectCores()

# the (negative) likelihood function which is minimized
# administrative right censoring
# 3 states
# based on Andersson & Keiding 2002
# States 1 and 2: Gompertz, State 3: Weibull - scale of Weibull is called rate here

loglikelihood <- function(data){
  state1 <- filter(data, to==1)
  state2 <- filter(data, to==2)
  state3 <- filter(data, to==3)
  if(length(state3$time)==0){
    f <- function(theta){
      shape <- theta[1:3]
      rate  <- theta[4:6]
      shape1 <- shape[1]
      shape2 <- shape[2]
      shape3 <- shape[3]
      rate1  <- rate[1]
      rate2  <- rate[2]
      rate3  <- rate[3]
      S <- function(t){exp(-(hgompertz(t,shape1,rate1)/shape1+hgompertz(t,shape2,rate2)/shape2+hweibull(t,shape3,rate3)*t/shape3)+rate1/shape1+rate2/shape2)}
      -(sum(log(S(data$time)))+sum(log(hgompertz(state1$time,shape1,rate1)))+sum(log(hgompertz(state2$time,shape2,rate2))))
    }
    return(f)
  }else{ 
    f <- function(theta){
      shape <- theta[1:3]
      rate  <- theta[4:6]
      shape1 <- shape[1]
      shape2 <- shape[2]
      shape3 <- shape[3]
      rate1  <- rate[1]
      rate2  <- rate[2]
      rate3  <- rate[3]
      S <- function(t){exp(-(hgompertz(t,shape1,rate1)/shape1+hgompertz(t,shape2,rate2)/shape2+hweibull(t,shape3,rate3)*t/shape3)+rate1/shape1+rate2/shape2)}
      -(sum(log(S(data$time)))+sum(log(hgompertz(state1$time,shape1,rate1)))+sum(log(hgompertz(state2$time,shape2,rate2)))+sum(log(hweibull(state3$time,shape3,rate3))))
    }
    return(f)
  }
}

# dat1, dat2 are the two datasets (for the two groups)
# joint likelihood is just the sum of the likelihoods 

joint_loglikelihood <- function(theta){
  theta1 <- theta[1:6]
  theta2 <- theta[7:12]
  loglikelihood(dat1)(theta1)+loglikelihood(dat2)(theta2)
}

# the constraint for the optimization procedure explained in the Algorithm
const <- function(theta){
  theta1 <- theta[1:6]
  theta2 <- theta[7:12]
  max(max(abs(hgompertz(grid,theta1[1],theta1[4])-hgompertz(grid,theta2[1],theta2[4]))),max(abs(hgompertz(grid,theta1[2],theta1[5])-hgompertz(grid,theta2[2],theta2[5]))),max(abs(hweibull(grid,theta1[3],theta1[6])-hweibull(grid,theta2[3],theta2[6]))))-delta
}

# underlying parameters of the study (based on the application examnple)
# corresponds to Scenario 3 in the paper (Table 3)
# for Scenario 4 just choose shape.dat2=shape.dat1 and rate.dat2=rate.dat1

# group 1
shape.gomp.dat<-c(-0.016359,-0.03582) 
rate.gomp.dat<-c(0.001906,0.00326) 
shape.weib.dat<-c(1.102)
scale.weib.dat<-c(2894.757)
shape.dat1 <- c(shape.gomp.dat,shape.weib.dat)
rate.dat1 <- c(rate.gomp.dat,scale.weib.dat)

# group 2
shape.gomp.dat2 <- c(-0.018423,-0.04259) 
rate.gomp.dat2 <- c(0.001606,0.00607) 
shape.weib.dat2<-c(1.114)
scale.weib.dat2<-c(1242.048)
shape.dat2 <- c(shape.gomp.dat2,shape.weib.dat2)
rate.dat2 <- c(rate.gomp.dat2,scale.weib.dat2)

grid <- seq(0.01,90,0.05) #for searching the maximum

t1 <- max(abs(hgompertz(grid,shape.dat1[1],rate.dat1[1])-hgompertz(grid,shape.dat2[1],rate.dat2[1])))
t2 <- max(abs(hgompertz(grid,shape.dat1[2],rate.dat1[2])-hgompertz(grid,shape.dat2[2],rate.dat2[2])))
t3 <- max(abs(hweibull(grid,shape.dat1[3],rate.dat1[3])-hweibull(grid,shape.dat2[3],rate.dat2[3])))
t <- max(t1,t2,t3)

deltav <- c(0.0028,0.004,0.005,0.007,0.01) 

# for saving results
rej <- matrix(NA,nrow=Nsim,ncol=length(deltav))
boot <- vector()

#### we are now ready to simulate

n.cores <- min(35, parallel::detectCores() - 1)
my.cluster <- parallel::makeCluster(n.cores, outfile="")

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

for(i in 1:length(deltav)){
  delta <- deltav[i]
  result <- foreach(j=1:Nsim,.combine=c,.packages=c("flexsurv","alabama","dplyr")) %dopar% {  
    set.seed(12345+j)
    dat1 <- sim.csh.gompertzandweibull(n=n1,
                                       shape.gomp=shape.gomp.dat,
                                       rate.gomp=rate.gomp.dat,
                                       shape.weib = shape.weib.dat,
                                       scale.weib = scale.weib.dat,
                                       cens.time=90)
    
    dat2 <- sim.csh.gompertzandweibull(n=n2,
                                       shape.gomp=shape.gomp.dat2,
                                       rate.gomp=rate.gomp.dat2,
                                       shape.weib = shape.weib.dat2,
                                       scale.weib = scale.weib.dat2,
                                       cens.time=90)
    
    #we estimate the transition intensities
    theta1 <- optim(par=c(shape.dat1,rate.dat1),fn=loglikelihood(dat1))$par
    theta2 <- optim(par=c(shape.dat2,rate.dat2),fn=loglikelihood(dat2))$par
    
    # calculating the test statistics
    t.stat1 <- max(abs(hgompertz(grid,theta1[1],theta1[4])-hgompertz(grid,theta2[1],theta2[4])))
    t.stat2 <- max(abs(hgompertz(grid,theta1[2],theta1[5])-hgompertz(grid,theta2[2],theta2[5])))
    t.stat3 <- max(abs(hweibull(grid,theta1[3],theta1[6])-hweibull(grid,theta2[3],theta2[6])))
    t.stat <- max(t.stat1,t.stat2,t.stat3)
    
    ### Constrained optimization  ###
    
    if (t.stat>=delta){
      min_cons <- c(theta1,theta2)}else{
        min_cons <- auglag(par=c(theta1,theta2),fn=joint_loglikelihood,heq=const,control.outer=list(method="nlminb",itmax=6,trace=FALSE))$par
      }
    theta1_cons <- min_cons[1:6]
    theta2_cons <- min_cons[7:12]
    
    # bootstrap (create new datasets with the estimated parameters)

    for(k in 1:B){
      dat1b <- sim.csh.gompertzandweibull(n=n1,
                                          shape.gomp=theta1_cons[1:2],
                                          rate.gomp=theta1_cons[4:5],
                                          shape.weib = theta1_cons[3],
                                          scale.weib = theta1_cons[6],
                                          cens.time=90)
      
      dat2b <- sim.csh.gompertzandweibull(n=n2,
                                          shape.gomp=theta2_cons[1:2],
                                          rate.gomp=theta2_cons[4:5],
                                          shape.weib = theta2_cons[3],
                                          scale.weib = theta2_cons[6],
                                          cens.time=90)
      theta1b <- NA
      theta2b <- NA
      try(theta1b <- optim(par=theta1_cons,fn=loglikelihood(dat1b))$par)
      try(theta2b <- optim(par=theta2_cons,fn=loglikelihood(dat2b))$par)
      
      if(any(is.na(theta1b)) | any(is.na(theta2b))){next}
      
      t.stat1b <- max(abs(hgompertz(grid,theta1b[1],theta1b[4])-hgompertz(grid,theta2b[1],theta2b[4])))
      t.stat2b <- max(abs(hgompertz(grid,theta1b[2],theta1b[5])-hgompertz(grid,theta2b[2],theta2b[5])))
      t.stat3b <- max(abs(hweibull(grid,theta1b[3],theta1b[6])-hweibull(grid,theta2b[3],theta2b[6])))
      boot[k] <- max(t.stat1b,t.stat2b,t.stat3b)  
    }
    
    # decision rule
    crit_val <- quantile(boot,alpha,na.rm=TRUE)
    #pval <- ecdf(boot)(t.stat1)
    if(t.stat<crit_val){
      rej=1
    }else{rej=0}
    
    print(j)
    return(rej)
  }
  rej[,i] <- result
  print(sum(result)/Nsim)
}

# evaluating the results
print(colSums(rej,na.rm=TRUE)/(Nsim-colSums(is.na(rej))))
