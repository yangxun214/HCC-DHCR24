library(survival)
library(survminer)
library(edgeR)
library(TCGAbiolinks)
library(dplyr)
head(myeloma)


load("DHCR24-TCGA-suvival.RData")
myeloma <-meta[,c(1:4)]

res.cut <- surv_cutpoint(myeloma, 
                         time = "time",
                         event = "event",
                         variables = c("DHCR24") 
)

summary(res.cut)
plot(res.cut, "DHCR24", palette = "npg")

res.cat <- surv_categorize(res.cut)
head(res.cat)
NUMBLE <- table(res.cat$DHCR24)
print(NUMBLE)

fit <- survfit(Surv(time, event) ~DHCR24, data = res.cat)#拟合生存分析

ggsurvplot(fit, 
           xlab="Months" ,
           ylab="Overall survival" ,
           pval=TRUE,
           pval.coord = c(5, 0.3),
           pval.size =5,
           conf.int=F, 
           palette = c("red", "black"),
           title="TCGA-LIHC", 
           legend.labs=c("DHCR24 High(n=170)", "DHCR24 Low(n=201)"), 
           legend.title=" ", 
           legend=c(0.3, 0.2)
)

