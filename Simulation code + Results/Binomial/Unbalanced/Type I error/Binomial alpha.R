

unlink(".RData")
unlink(".Rhistory")


library(multcomp)
library(lme4)
library(Rcpp)
library(arm)
library(brglm)


Simfun_all <- function(n_center1, n_center2, n_center3, n_center4,
                       n_center5, p_center1, p_center2,
                       p_center3, p_center4, p_center5, nsim, norep=norep){
  
  
  coverage_BayesGLM <- logical(length=nsim)
  coverage_BrGLM <- logical(length=nsim)
  coverage_GLM <- logical(length=nsim)
  zeros <- logical(length=nsim)
  type1_GLM <- logical(length=nsim)
  type1_BayesGLM <- logical(length=nsim)
  type1_BrGLM <- logical(length=nsim)
  
  
   nc1 <- round(n_center1/2, digits = 0)
  nc2 <- round(n_center2/5, digits = 0)
  nc3 <- round(n_center3, digits = 0)
  nc4 <- round(n_center4/4, digits = 0)
  nc5 <- round(n_center5*2, digits = 0)
  
  
  for(i in 1:nsim){
    
    
    center1 <- rbinom(n=norep, size =nc1 , prob=p_center1)
    center2 <- rbinom(n=norep, size =nc2 , prob=p_center2)
    center3 <- rbinom(n=norep, size =nc3 , prob=p_center3)
    center4 <- rbinom(n=norep, size =nc4 , prob=p_center4)
    center5 <- rbinom(n=norep, size=nc5,  prob=p_center5)
    
    succ <- c(center1, center2, center3, center4, center5)
    size <- rep(c(nc1, nc2, nc3, nc4, nc5), each=norep)
  #  print(size)
    
    
    label <- rep(c("center1", "center2", "center3", "center4", "center5"), each=norep)
    center <- as.factor(label)
    fail <- size - succ
    
    n <- rep(c(nc1, nc2, nc3, nc4, nc5), each=norep)
    p <- rep(c(p_center1, p_center2, p_center3, p_center4, p_center5), each=norep)
    
    
    dat <- data.frame(center=center, n=n, p=p, succ=succ, fail=fail)
     print(dat)
    
    fit_BayesGLM <- bayesglm(cbind(succ, fail) ~ center, data=dat, family=binomial)
    compGrandMean_BayesGLM <- glht(fit_BayesGLM, mcp(center="GrandMean"))
    center5f_interval_BayesGLM <- confint(compGrandMean_BayesGLM)
    confintervals_BayesGLM <- as.data.frame(center5f_interval_BayesGLM$confint)
    confintervals_BayesGLM$method <- c("BayesGLM")
   # print(confintervals_BayesGLM)
    
    scompGrandMean_BayesGLM <- summary(compGrandMean_BayesGLM)
    pvalues_BayesGLM <- scompGrandMean_BayesGLM$test$pvalues
    compnames_BayesGLM <- names(scompGrandMean_BayesGLM$test$coefficients)
    
    tstres_BayesGLM <- data.frame(centerment=compnames_BayesGLM,
                                  pvalues_BayesGLM=pvalues_BayesGLM)
    
    
    
    fit_BrGLMGLM <- brglm(cbind(succ, fail) ~ center, data=dat, family=binomial)
    compGrandMean_BrGLM <- glht(fit_BrGLMGLM, mcp(center="GrandMean"))
    center5f_interval_BrGLM <- confint(compGrandMean_BrGLM)
    confintervals_BrGLM <- as.data.frame(center5f_interval_BrGLM$confint)
    confintervals_BrGLM$method <- c("BrGLM")
  #  print(confintervals_BrGLM)
    scompGrandMean_BrGLM <- summary(compGrandMean_BrGLM)
    pvalues_BrGLM <- scompGrandMean_BrGLM$test$pvalues
    compnames_BrGLM <- names(scompGrandMean_BrGLM$test$coefficients)
    tstres_BrGLM <- data.frame(centerment=compnames_BrGLM,
                               pvalues_BrGLM=pvalues_BrGLM)
    
    
    
    
    fitGLM <- glm(cbind(succ, fail) ~ center, data=dat, family=binomial)
    compGrandMeanGLM <- glht(fitGLM, mcp(center="GrandMean"))
    center5f_interval_GLM <- confint(compGrandMeanGLM)
    confintervals_GLM <- as.data.frame(center5f_interval_GLM$confint)
    confintervals_GLM$method <- c("GLM")
    
    scompGrandMean_GLM <- summary(compGrandMeanGLM)
    pvalues_GLM <- scompGrandMean_GLM$test$pvalues
   #  print(pvalues_GLM)
    compnames_GLM <- names(scompGrandMean_GLM$test$coefficients)
    
    tstres_GLM <- data.frame(centerment=compnames_GLM,
                             pvalues_GLM=pvalues_GLM)
    
    
    # print(tstres_GLM$pvalues)
    #print(confintervals_GLM)
    
    logit <- function(p){log(p)-log(1-p)}
    ni <- c(n_center5, n_center1, n_center2, n_center3, n_center4)
    pi <- c(p_center5, p_center1, p_center2, p_center3, p_center4)
    cm <- contrMat(n=ni, type="GrandMean")
   # print(cm)
    logitpi <- matrix(logit(pi), ncol=1)
    truelogoddsratio <- cm %*% logitpi
   # print(truelogoddsratio)
    
    coverage_BayesGLM_CI <- all({confintervals_BayesGLM$lwr[confintervals_BayesGLM$method=="BayesGLM"]<= truelogoddsratio & confintervals_BayesGLM$upr[confintervals_BayesGLM$method=="BayesGLM"]>= truelogoddsratio})
    coverage_BayesGLM[i] <- coverage_BayesGLM_CI
   # print(coverage_BayesGLM)
    coverage_BrGLM_CI <- all({confintervals_BrGLM$lwr[confintervals_BrGLM$method=="BrGLM"]<= truelogoddsratio & confintervals_BrGLM$upr[confintervals_BrGLM$method=="BrGLM"]>= truelogoddsratio})
    coverage_BrGLM[i] <- coverage_BrGLM_CI
   # print(coverage_BrGLM)
    coverage_GLM_CI <- all({confintervals_GLM$lwr[confintervals_GLM$method=="GLM"]<= truelogoddsratio & confintervals_GLM$upr[confintervals_GLM$method=="GLM"]>= truelogoddsratio})
    coverage_GLM[i] <-   coverage_GLM_CI
   # print(coverage_GLM)
    
    
    
    rej_GLM <- any(tstres_GLM$pvalues_GLM < 0.05)
    #print(tstres_GLM$pvalues_GLM)
    #print(rej_GLM)
    type1_GLM[i] <- rej_GLM
    # print(type1_GLM)
    
    rej_BrGLM <- any(tstres_BrGLM$pvalues_BrGLM < 0.05)
    # print(tstres_BrGLM$pvalues_BrGLM )
    type1_BrGLM[i] <- rej_BrGLM
    
    
    rej_BayesGLM <- any(tstres_BayesGLM$pvalues_BayesGLM < 0.05)
    # print(tstres_BayesGLM$pvalues_BayesGLM)
    type1_BayesGLM[i] <- rej_BayesGLM
    
    
    
    zerodat <- subset(dat, center=="center5")
    #print(zerodat)
    #print(zerodat)
    zeros[i] <- all({zerodat$succ==0})
    #print(zeros)
    
  }
  est.coverageprob_BayesGLM <- sum(coverage_BayesGLM)/nsim
  est.alpha_BayesGLM <- sum(type1_BayesGLM, na.rm=TRUE)/nsim
  prop.fail_BayesGLM <- sum(is.na(type1_BayesGLM))/nsim
  
  
  
  est.coverageprob_BrGLM <- sum(coverage_BrGLM)/nsim
  est.alpha_BrGLM <- sum(type1_BrGLM, na.rm=TRUE)/nsim
  prop.fail_BrGLM <- sum(is.na(type1_BrGLM))/nsim
  
  
  est.coverageprob_GLM <- sum(coverage_GLM)/nsim
  est.alpha_GLM <- sum(type1_GLM, na.rm=TRUE)/nsim
  prop.fail_GLM <- sum(is.na(type1_GLM))/nsim
  
  
  
  
  zerodatcov <- sum(zeros)/nsim
  
  
  result <- data.frame(n_center5=n_center5,
                       n_center1=n_center1,
                       n_center2=n_center2,
                       n_center3=n_center3,
                       n_center4=n_center4,
                       p_center5=p_center5, 
                       p_center1=p_center1,
                       p_center2=p_center2,
                       p_center3=p_center3,
                       p_center4=p_center4,
                       norep=norep,
                       CI_BayesGLM= est.coverageprob_BayesGLM,
                       CI_BrGLM=est.coverageprob_BrGLM,
                       CI_GLM= est.coverageprob_GLM,
                       est.alpha_BayesGLM=est.alpha_BayesGLM,
                       est.alpha_BrGLM=est.alpha_BrGLM,
                       est.alpha_GLM=est.alpha_GLM,
                       zerodatcov=zerodatcov,
                       prop.fail_BayesGLM=prop.fail_BayesGLM,
                       prop.fail_BrGLM=prop.fail_BrGLM,
                       prop.fail_GLM=prop.fail_GLM,
                       nsim=nsim)
  
  return(result)
}

# Simfun_all(n_center1=5, n_center2=5, n_center3=5, n_center4=5,
# n_center5=5, p_center1=0.2, p_center2=0.2,p_center3=0.2, p_center4=0.2, p_center5=0.2, nsim=1, norep=1)



simdat <- data.frame("n"=c(5,10,
                           20, 25,30,
                           50,80,
                           100, 200,
                           150,300,
                           400,
                           500, 1000,
                           1500, 2000), "E1"=c(2, 2,
                                               5,9,
                                               23,44,
                                               3,33,
                                               24,5,
                                               45,11,
                                               89, 180,
                                               200, 280),
                     "E2"=c(2, 2,
                            2,3,
                            3,3,
                            3,4,
                            5,5,
                            10,11,
                            22, 18,
                            53, 80),
                     "j"=1)





# simdat <- data.frame("n"=c(50), "E1"=c(8),
#                      "E2"=c(29),
# 
#                      "E3"=c(44),
#                      "j"=1)
#simdat$p3 <- simdat$E3/(simdat$n*simdat$j)

simdat$p1 <- simdat$E1/(simdat$n*simdat$j)
  simdat$p2 <- simdat$E2/(simdat$n*simdat$j)
  



system.time(
  
  datfra <- apply(simdat, MARGIN = 1, 
                  function(x){
                    Simfun_all(nsim=1000,
                               n_center1=x["n"],
                               n_center2=x["n"],
                               n_center3=x["n"],
                               n_center4=x["n"],
                               n_center5=x["n"],
                               p_center1=x["p1"],
                               p_center2=x["p1"],
                               p_center3=x["p1"],
                               p_center4=x["p1"],
                               p_center5=x["p1"], 
                               norep=x["j"]
                    )})
  
)


fsimdat <-do.call(rbind.data.frame, datfra)
fsimdat[,-c(1:4)]

write.csv2(fsimdat, "Binomial unbalanced alpha_1000.csv")







