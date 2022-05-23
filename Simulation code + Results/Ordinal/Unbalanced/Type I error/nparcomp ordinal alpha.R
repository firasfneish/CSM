

unlink(".RData")
unlink(".Rhistory")

library(ANOM)
library(MCPAN)
library(nparcomp)
library(tidyverse)




datgen <- function(n1, n2, mean1, sd1, mean2, sd2){

  nc1 <- n1#sample(n1:n2, 1)
  nc2 <- n1#sample(n1:n2, 1)
  nc3 <- n1#sample(n1:n2, 1)
  nc4 <- n1#sample(n1:n2, 1)
  nc5 <- n1#sample(n1:n2, 1)
  nc6 <- n1#sample(n1:n2, 1)
  nc7 <- n1#sample(n1:n2, 1)
  nc8 <- n1#sample(n1:n2, 1)
  nc9 <- n1#sample(n1:n2, 1)
  nc10 <- n2#sample(n3:n4, 1)
  

    c1 <- round(rnorm(n=nc1, mean=mean1, sd=sd1), digits = 0)
  c2 <-  round(rnorm(n=nc2, mean=mean1, sd=sd1), digits = 0)
  c3 <-  round(rnorm(n=nc3,  mean=mean1, sd=sd1), digits = 0)
  c4 <- round(rnorm(n=nc4, mean=mean1, sd=sd1), digits = 0)
  c5 <- round( rnorm(n=nc5, mean=mean1, sd=sd1), digits = 0)
  c6 <- round( rnorm(n=nc6,  mean=mean1, sd=sd1), digits = 0)
  c7 <- round( rnorm(n=nc7, mean=mean1, sd=sd1), digits = 0)
  c8 <- round( rnorm(n=nc8, mean=mean1, sd=sd1), digits = 0)
  c9 <- round( rnorm(n=nc9, mean=mean1, sd=sd1), digits = 0)
  c10 <- round( rnorm(n=nc10, mean=mean2, sd=sd2), digits = 0)

  dat <- qpcR:::cbind.na(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10) %>% as.data.frame()
#print(dat)
  dat <- dat %>% pivot_longer(names_to = "centers", values_to="val", 1:10) %>%
    mutate(centers=as.factor(centers)) %>% filter(!is.na(val))


  obser_per_centre <- dat %>% group_by(centers) %>% summarise(n=n())
  dat1 <- dat %>% left_join(obser_per_centre, by="centers")

  return(dat1)
}


  

pvalgen <- function(dat){
ss <- tapply(dat$val, dat$centers, length)
#print(ss)
Mat <- contrMat(ss, type="GrandMean")
wf3 <- mctp(val ~ centers, data=dat, type="UserDefined",
            contrast.matrix=Mat, alternative="two.sided", info=FALSE,
            correlation=TRUE, asy.method="mult.t")
#print(wf3)
#ANOM(wf3)

wf4 <- wf3$Analysis %>% rownames_to_column("Centers") %>% cbind( wf3$Data.Info$Size)
return(wf4)
}




simfun <- function(n1, n2,  mean1, sd1, mean2, sd2, nsim){

  type1 <- logical(length=nsim)

  for(i in 1:nsim){

    dat <- datgen(n1=n1, n2=n2, mean1 = mean1, mean2=mean2, sd1=sd1, sd2=sd2)
#print(dat)
    pvals <- pvalgen(dat)
    nc1 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 1"]
    ## this is not a mistake, its the naming order only
    nc10 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 2"]

    nc2 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 3"]
    nc3 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 4"]
    nc4 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 5"]
    nc5 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 6"]
    nc6 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 7"]
    nc7 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 8"]
    nc8 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 9"]
    nc9 <- pvals$`wf3$Data.Info$Size`[pvals$Centers=="C 10"]

    rej <- any(pvals$p.Value < 0.05)
    #print(rej)

    type1[i] <- rej

      }


    est.alpha <- sum(type1)/nsim
  #print(est.alpha)
    pval <- data.frame(n1=n1,
                     n2=n2,
                     mean1=mean1,
                     mean2=mean2,
                     sd1=sd1,
                     sd2=sd2,
                     est.alpha=est.alpha,
                     nsim=nsim)

  return(pval)
  }


#simfun(n1=10, n2=3, mean1 = 35, mean2 = 35, sd1=5, sd2=0,  nsim = 1)



# simdat <- expand.grid(n1=c(3, 5, 10, 20, 40, 50, 80),
#                       n2=c(15,25,100, 200),
#                       n3=c(2, 4),
#                       n4=c(3, 7, 10), mean1 = 35, mean2 = 35, sd1=5, sd2=0)



nm= c(3, 5, 10, 20, 40, 50, 80, 100, 200)



n1=nm*2
n2=nm

simdat <- data.frame(n1=n1, 
                     n2=n2, 
                     
                     mean1=35,  mean2=35,   sd1=1, sd2=1) 


system.time(
  datfradiff <- apply(simdat, MARGIN = 1,
                      function(x){
                        simfun(nsim=1000,
                               n1=x["n1"],
                               n2=x["n2"],
                               mean1=x["mean1"],
                               mean2=x["mean2"],
                               sd1=x["sd1"],
                               sd2=x["sd2"]
                              )})
  )

datfradiff

da <- data.frame(datfradiff)
da

fsimdat <-do.call(rbind.data.frame, datfradiff) 

write.csv2(fsimdat, "nparcom ordinal unbalanced alpha_1000.csv")



