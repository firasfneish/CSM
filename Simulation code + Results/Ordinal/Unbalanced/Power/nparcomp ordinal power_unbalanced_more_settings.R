

unlink(".RData")
unlink(".Rhistory")

library(ANOM)
library(MCPAN)
library(nparcomp)
library(tidyverse)
library(pwr)



datgen <- function(n1, n2,  mean1, d, sd1, sd2){

  mean2= mean1-d
  
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




simfun <- function(n1, n2,  mean1, sd1, d, sd2, nsim){

  type1 <- logical(length=nsim)

  for(i in 1:nsim){

    dat <- datgen(n1=n1, n2=n2,  mean1 = mean1, d=d, sd1=sd1, sd2=sd2)
#print(dat)
    pvals <- pvalgen(dat)
    

    rej <- any(pvals$p.Value < 0.05)
    #print(rej)

    type1[i] <- rej

      }


    est.alpha <- sum(type1)/nsim
  #print(est.alpha)
    pval <- data.frame(n1=n1,
                     n2=n2,
                     mean1=mean1,
                     d=d,
                     sd1=sd1,
                     sd2=sd2,
                     est.alpha=est.alpha,
                     nsim=nsim)

  return(pval)
  }



nm= c(2, 3, 5, 10, 20, 40, 50, 80, 100, 200)

n1=nm*2
n2=nm


simdat1 <- data.frame(n1=n1, 
                     n2=n2, 
                     
                     mean1=35,   sd1=1, sd2=1) %>%   mutate(
                       d=c(pwr.t2n.test(n1 =n1[1]  , n2= n2[1], sig.level =0.05, power = 0.8)$d,
                           pwr.t2n.test(n1 =n1[2]  , n2= n2[2], sig.level =0.05, power = 0.8)$d, 
                           pwr.t2n.test(n1 =n1[3]  , n2= n2[3], sig.level =0.05, power = 0.8)$d,
                           pwr.t2n.test(n1 =n1[4]  , n2= n2[4], sig.level =0.05, power = 0.8)$d, 
                           pwr.t2n.test(n1 =n1[5]  , n2= n2[5], sig.level =0.05, power = 0.8)$d, 
                           pwr.t2n.test(n1 =n1[6]  , n2= n2[6], sig.level =0.05, power = 0.8)$d,
                           pwr.t2n.test(n1 =n1[7]  , n2= n2[7], sig.level =0.05, power = 0.8)$d,
                           pwr.t2n.test(n1 =n1[8]  , n2= n2[8], sig.level =0.05, power = 0.8)$d,
                           pwr.t2n.test(n1 =n1[9]  , n2= n2[9], sig.level =0.05, power = 0.8)$d, 
                           pwr.t2n.test(n1 =n1[10]  , n2= n2[10], sig.level =0.05, power = 0.8)$d) ) 


n1=nm*10
n2=nm

simdat2 <- data.frame(n1=n1, 
                      n2=n2, 
                      
                      mean1=35,   sd1=1, sd2=1) %>%   mutate(
                        d=c(pwr.t2n.test(n1 =n1[1]  , n2= n2[1], sig.level =0.05, power = 0.8)$d,
                            pwr.t2n.test(n1 =n1[2]  , n2= n2[2], sig.level =0.05, power = 0.8)$d, 
                            pwr.t2n.test(n1 =n1[3]  , n2= n2[3], sig.level =0.05, power = 0.8)$d,
                            pwr.t2n.test(n1 =n1[4]  , n2= n2[4], sig.level =0.05, power = 0.8)$d, 
                            pwr.t2n.test(n1 =n1[5]  , n2= n2[5], sig.level =0.05, power = 0.8)$d, 
                            pwr.t2n.test(n1 =n1[6]  , n2= n2[6], sig.level =0.05, power = 0.8)$d,
                            pwr.t2n.test(n1 =n1[7]  , n2= n2[7], sig.level =0.05, power = 0.8)$d,
                            pwr.t2n.test(n1 =n1[8]  , n2= n2[8], sig.level =0.05, power = 0.8)$d,
                            pwr.t2n.test(n1 =n1[9]  , n2= n2[9], sig.level =0.05, power = 0.8)$d, 
                            pwr.t2n.test(n1 =n1[10]  , n2= n2[10], sig.level =0.05, power = 0.8)$d) ) 

simdat <- simdat1 |> rbind(simdat2)

system.time(
  datfradiff <- apply(simdat, MARGIN = 1,
                      function(x){
                        simfun(nsim=1000,
                               n1=x["n1"],
                               n2=x["n2"],
                               mean1=x["mean1"],
                               d=x["d"],
                               sd1=x["sd1"],
                               sd2=x["sd2"]
                              )})
  )

datfradiff

da <- data.frame(datfradiff)
da

fsimdat <-do.call(rbind.data.frame, datfradiff) 
fsimdat


write.csv2(fsimdat, "nparcom ordinal unbalanced power_1000.csv")



