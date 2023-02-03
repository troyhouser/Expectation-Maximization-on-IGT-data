library(mclust)

setwd("/Users/troyhouser/Downloads/IGTdataSteingroever2014")
load("IGTdata.rdata")
X = cbind(lo_95[1,],wi_95[1,])
rownames(X) = NULL
fit = Mclust(X,G=4)
fit$classification
summary(fit,parameters=T)
plot(fit, what="BIC")
dens = densityMclust(X)
summary(dens)
dr = MclustDR(fit)
summary(dr)
plot(dr,what="boundaries",ngrid=200)
plot(X,col=predict(fit)$classification)
######expectation step
expectation = function(sample,p,a,b){
  p_expectation = (p*dbinom(sample,1,a)) / (p*dbinom(sample,1,a)+(1-p)*
                                              dbinom(sample,1,b))
  return(p_expectation)
}

#######maximization step
maximization = function(sample,epart){
  #estimate p
  p_temp = mean(epart)
  
  #estimate a and b
  a_temp = sum(sample*epart) / sum(epart)
  b_temp = sum(sample*(1-epart)) / sum(1-epart)
  
  list(p_temp,a_temp,b_temp)
}

########3 EM algorithm
EM = function(sample,p_inits,a_inits,b_inits,maxit=1000,tol=1e-6){
  ###estimation of initial parameters
  flag = 0
  p_cur = p_inits; a_cur = a_inits; b_cur = b_inits
  
  # iterate between expectation and maximization parts
  for(i in 1:maxit){
    cur = c(p_cur,a_cur,b_cur)
    new = maximization(sample,expectation(sample,p_cur,a_cur,b_cur))
    p_new = new[[1]]; a_new = new[[2]]; b_new = new[[3]]
    new_step = c(p_new,a_new,b_new)
    
    # stop iteration if difference between current and new estimates is less than tol
    if(all(abs(cur - new_step) < tol)){flag = 1; break}
    
    ## otherwise:
    p_cur = p_new; a_cur = a_new; b_cur = b_new
  }
  if(!flag) warning("didn't converge\n")
  
  list(p_cur,a_cur,b_cur)
}

infomat = function(sample,p.est,a.est,b.est){
  expectation.est = expectation(sample,p.est,a.est,b.est)
  info.mat = matrix(rep(0,9),3,3)
  info.mat[1,1] = sum(expectation.est) / (p.est^2) - 
    sum((1-expectation.est)) / ((1-p.est)^2)
  info.mat[2,2] = sum(expectation.est * sample) / (a.est^2) -
    sum(expectation.est * (1-sample)) / ((1-a.est)^2)
  info.mat[3,3] = sum((1-expectation.est) * sample) / (b.est^2) -
    sum((1-expectation.est) * (1-sample)) / ((1-b.est)^2)
}
freqs = as.vector(table(choice_95[1,]))
p_true = (freqs[1]+freqs[2]) / sum(freqs)
a_true = freqs[1] / sum(freqs[1],freqs[2])
b_true = freqs[2] / sum(freqs[1],freqs[2])
n <- 10000
true.out <- c(p_true,a_true,b_true)
u <- ifelse(runif(n)<p_true, rbinom(n,1,a_true),rbinom(n,1,b_true))

# Set parameter estimates
p_init = 0.50; a_init = 0.5; b_init = 0.50

output <- EM(u,p_init,a_init,b_init)

# Confidence Intervals
sd.out <- sqrt(diag(solve(infomat(u,output[[1]],output[[2]],output[[3]]))))
data.frame("Truth" = true.out, "EM Estimate" = unlist(output), "Lower CI" = unlist(output) - qnorm(.975)*sd.out, "Upper CI" = unlist(output) + qnorm(.975)*sd.out)


##########################################################################################

### multinomial hmm
c1 = c2 = c3 = c4 = 0
for(i in 1:ncol(choice_95)){
  if(choice_95[1,i]==1){
    c1 = c1 + 1
  }else if(choice_95[1,i]==2){
    c2 = c2 + 1
  }else if(choice_95[1,i]==3){
    c3 = c3 + 1
  }else if(choice_95[1,i]==4){
    c4 = c4 + 1
  }
}
p1 = c1/95
p2 = c2/95
p3 = c3/95
p4 = c4/95
my_prob = c(p1,p2,p3,p4)
experiments = rmultinom(n=nrow(choice_95),size=ncol(choice_95),prob=my_prob)
