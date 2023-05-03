

unlink(".RData")
unlink(".Rhistory")

library(ANOM)
library(MCPAN)
library(nparcomp)
library(tidyverse)
library(pwr)



Packages <- c("tidyverse", "multcomp", "arm", "brglm", "kableExtra", "ggpubr", "RColorBrewer")
lapply(Packages, library, character.only = TRUE)


datgen <- function(n1, n2,  mean1, sd1, d, sd2){

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

    c1 <- rnorm(n=n1, mean=mean1, sd=sd1)
  c2 <-  rnorm(n=n1, mean=mean1, sd=sd1)
  c3 <-  rnorm(n=n1,  mean=mean1, sd=sd1)
  c4 <- rnorm(n=n1, mean=mean1, sd=sd1)
  c5 <-  rnorm(n=n1, mean=mean1, sd=sd1)
  c6 <-  rnorm(n=n1,  mean=mean1, sd=sd1)
  c7 <-  rnorm(n=n1, mean=mean1, sd=sd1)
  c8 <-  rnorm(n=n1, mean=mean1, sd=sd1)
  c9 <-  rnorm(n=n1, mean=mean1, sd=sd1)
  c10 <-  rnorm(n=n2, mean=mean2, sd=sd2)

  dat <- qpcR:::cbind.na(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10) %>% as.data.frame()
#print(dat)
  dat <- dat %>% pivot_longer(names_to = "centers", values_to="val", 1:10) %>%
    mutate(centers=as.factor(centers)) %>% filter(!is.na(val))


  obser_per_centre <- dat %>% group_by(centers) %>% summarise(n=n())
  dat1 <- dat %>% left_join(obser_per_centre, by="centers")

  return(dat1)
}

#dat <- datgen(n1=5, n2=90,  mean1=35, sd1=1, d=2, sd2=1)

pvalgen <- function(dat){
  fit <- lm(val ~ centers, data = dat)
  compGrandMean <- glht(fit, mcp(centers="GrandMean"))
  center5f_interval <- confint(compGrandMean)
  confintervals <- as.data.frame(center5f_interval$confint)
  scompGrandMean <- summary(compGrandMean)
  p.values <- scompGrandMean$test$pvalues
 lm_out <- confintervals %>% cbind(p.values)%>% rownames_to_column("Centers") 
 

 ss <- tapply(dat$val, dat$centers, length)

 Mat <- contrMat(ss, type="GrandMean")
 wf3 <- mctp(val ~ centers, data=dat, type="UserDefined",
             contrast.matrix=Mat, alternative="two.sided", info=FALSE,
             correlation=TRUE, asy.method="mult.t")


wf5 <- data.frame(centers=lm_out$Centers ,nparcomp_pval=wf3$Analysis$p.Value, lm_pval=lm_out$p.values )

 return(wf5)
 
 
}



#pvalgen(dat)


simfun <- function(n1, n2,  mean1, sd1, d, sd2, nsim){

  type1_lm <- logical(length=nsim)
  type1_nparcomp <- logical(length=nsim)
  
  for(i in 1:nsim){
    dat <- datgen(n1=n1, n2=n2, mean1 = mean1, d=d, sd1=sd1, sd2=sd2)
    pvals <- pvalgen(dat)
    rej_lm <- any(pvals$lm_pval < 0.05)
    rej_nparcomp <- any(pvals$nparcomp_pval < 0.05)
    type1_lm[i] <- rej_lm

    type1_nparcomp[i] <- rej_nparcomp
    
      }


    est.alpha_lm <- sum(type1_lm)/nsim

    est.alpha_nparcomp <- sum(type1_nparcomp)/nsim
    
    pval <- data.frame(n1=n1,
                     n2=n2,
                   
                     mean1=mean1,
                     d=d,
                     sd1=sd1,
                     sd2=sd2,
                     est.alpha_lm=est.alpha_lm,
                     est.alpha_nparcomp=est.alpha_nparcomp,
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
                         pwr.t2n.test(n1 =n1[9]  , n2= n2[9], sig.level =0.05, power = 0.8)$d  , 
                         pwr.t2n.test(n1 =n1[10]  , n2= n2[10], sig.level =0.05, power = 0.8)$d  ) )



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
                           pwr.t2n.test(n1 =n1[10]  , n2= n2[10], sig.level =0.05, power = 0.8)$d   ) )


simdat <-  simdat1 |> rbind(simdat2) 

system.time(
  datfradiff <- apply(simdat, MARGIN = 1000,
                      function(x){
                        simfun(nsim=1,
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

write.csv2(fsimdat, "lm_nparcomp unbalanced power_1000.csv")

nparcomp_unbalanced_power <- fsimdat
xaxislabel=10
yaxislabel=10
titleyaxis=5
titlexaxis=5
points_size=5
order_legend=5
legend.text=5
hline=3





gg_nparcomp_power_unbalanced <- nparcomp_unbalanced_power %>% pivot_longer(names_to = "Method", values_to="est.alpha", c(est.alpha_lm, est.alpha_nparcomp)) %>% mutate(Method= case_when(Method== "est.alpha_lm" ~ "Linear Model", Method=="est.alpha_nparcomp"~"Non-parametric")) %>%
  mutate(n1=factor(n1, levels = factor(unique(sort(nparcomp_unbalanced_power$n1))) ), 
         n2=factor(n2, levels = factor(unique(sort(nparcomp_unbalanced_power$n2))) )) %>% 
  
  ggplot(aes(n1, est.alpha))+theme_bw()+  
  geom_point(position=position_jitter(width=0,height=0), aes(shape=(n2), stroke=2, color=d ), size=points_size)+
  theme(axis.text.x = element_text(angle=0) ) +geom_hline(yintercept=0.8, linetype="dashed", 
                                                          color = "red", size=hline)+theme(strip.text = element_text(size=19), 
                                                                                           axis.text.x = element_text(size = xaxislabel),
                                                                                           axis.text.y = element_text(size = yaxislabel))+ #ggtitle("Unbalanced")+
  scale_shape_manual(values=c(1:16), name=expression(bold(n[i] ~ "in deviating center")))+
  labs(y = expression(1-beta))+
  theme(legend.title = element_text(size=legend.text+20+15+20, face = "bold"), legend.text = element_text(size=legend.text+50, face = "bold"), axis.title.y = element_text(size=30))+
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, size=30))+ geom_line(aes(group=Method, linetype=Method), size=1.5)+theme(axis.title.x = element_text(size=19), legend.position="bottom")+xlab(expression(bold(n[i] ~ "per center")))+theme(axis.title.y = element_text(size = titleyaxis), 
                                                                                                                                                                                                                                          axis.title.x = element_text(size = titlexaxis))+theme(legend.spacing = unit(5, unit="cm"))+
  scale_color_continuous(
    limits=c(0,3.1), breaks=c(0.2, 1.5, 2.5),
    guide = guide_colourbar(barwidth = 40, barheight = 8, order = 1,title.vjust = 1, label.theme = element_text(size = 65, face = "bold"), title.theme = element_text(size=150, face = "bold")), low="blue",  high="red",
    name=expression(delta)) + scale_linetype_manual(values=c("solid", "dotted"))+   guides(linetype=guide_legend(keywidth = 5, direction = "vertical") )



gg_nparcomp_power_unbalanced 
