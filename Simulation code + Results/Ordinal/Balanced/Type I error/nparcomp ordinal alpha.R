

unlink(".RData")
unlink(".Rhistory")

library(ANOM)
library(MCPAN)
library(nparcomp)
library(tidyverse)




datgen <- function(n, mean1, sd1, mean2, sd2){


    c1 <- round(rnorm(n=n, mean=mean1, sd=sd1), digits = 0)
  c2 <-  round(rnorm(n=n, mean=mean1, sd=sd1), digits = 0)
  c3 <-  round(rnorm(n=n,  mean=mean1, sd=sd1), digits = 0)
  c4 <- round(rnorm(n=n, mean=mean1, sd=sd1), digits = 0)
  c5 <- round( rnorm(n=n, mean=mean1, sd=sd1), digits = 0)
  c6 <- round( rnorm(n=n,  mean=mean1, sd=sd1), digits = 0)
  c7 <- round( rnorm(n=n, mean=mean1, sd=sd1), digits = 0)
  c8 <- round( rnorm(n=n, mean=mean1, sd=sd1), digits = 0)
  c9 <- round( rnorm(n=n, mean=mean1, sd=sd1), digits = 0)
  c10 <- round( rnorm(n=n, mean=mean2, sd=sd2), digits = 0)

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




simfun <- function(n, mean1, sd1, mean2, sd2, nsim){

  type1 <- logical(length=nsim)

  for(i in 1:nsim){

    dat <- datgen(n=n, mean1 = mean1, mean2=mean2, sd1=sd1, sd2=sd2)
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
    pval <- data.frame(n=n,
                     mean1=mean1,
                     mean2=mean2,
                     sd1=sd1,
                     sd2=sd2,
                     est.alpha=est.alpha,
                     nsim=nsim)

  return(pval)
  }


#simfun(n=10, mean1 = 35, mean2 = 35, sd1=5, sd2=0,  nsim = 1)


# 
# simdat <- expand.grid(n1=c(3, 5, 10, 20, 40, 50, 80),
#                       n2=c(15,25,100, 200),
#                       n3=c(2, 4),
#                       n4=c(3, 7, 10), mean1 = 35, mean2 = 35, sd1=5, sd2=0)



simdat <- expand.grid(n=c(3, 5, 10, 20, 40, 50, 80, 100, 200),
                      mean1 = 35, mean2 = 35, sd1=5, sd2=0)



system.time(
  datfradiff <- apply(simdat, MARGIN = 1,
                      function(x){
                        simfun(nsim=1000,
                               n=x["n"],
                               mean1=x["mean1"],
                               mean2=x["mean2"],
                               sd1=x["sd1"],
                               sd2=x["sd2"]
                              )})
  )

datfradiff

da <- data.frame(datfradiff)
da

fsimdat <-do.call(rbind.data.frame, datfradiff) #%>% #select(-c(nc1, nc2, nc3, nc4, nc5, nc6, nc7, nc8, nc9, nc10)) %>%
  #mutate(N_each=paste(n1,n2,n3, n4, sep = ":")) %>% mutate(group=case_when(est.alpha>0 ~ "Maybe problem", TRUE ~ "No problem"))
write.csv2(fsimdat, "nparcom ordinal alpha_1000_balanced.csv")

# fsimdat %>% filter(group %in% c("Maybe problem"))%>% ggplot(aes(N_each, est.alpha))+  geom_point(position=position_jitter(width=0.3,height=0))+
#   theme(axis.text.x = element_text(angle=90) ) +facet_wrap(~ group)


