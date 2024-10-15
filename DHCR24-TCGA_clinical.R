library(survival)
library(plyr)

load("DHCR24_TCGA-clinical.RData")
str(aa)


#=========Univariable Regression Model===============
aa$event<-factor(aa$event)
summary(aa$event)
y<- Surv(time=aa$time,event=aa$event==1)

Uni_cox_model<-function(x){
  FML<-as.formula(paste0("y~",x))
  Cox<-coxph(FML,data = aa) 
  Sum<-summary(Cox)
  CI<-paste0(round(Sum$conf.int[,3:4],3),collapse = "-") 
  Pvalue<-round(Sum$coefficients[,5],3)
  HR<-round(Sum$coefficients[,2],3)
  Unicox<-data.frame("Characteristics"=x,
                     "Hazard Ratio"=HR,
                     "CI95"=CI,
                     "P value"=Pvalue)
  return(Unicox)
}

names(aa)
variable.names<- colnames(aa)[c(5:16)] 

Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox<- ldply(Uni_cox,data.frame)
Uni_cox$CI<-paste(Uni_cox$CI5,'-',Uni_cox$CI95)
Uni_cox<-Uni_cox[,-5]


#=========Multiple Regression Model===============
Uni_cox$Characteristics[Uni_cox$P.value<0.05]

mul_cox_model<- as.formula(paste0 ("y~",
                                   paste0(Uni_cox$Characteristics[Uni_cox$P.value<0.05],
                                          collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=aa)
cox4<-summary(mul_cox) 


mul_HR<- round(cox4$coefficients[,2],3) 
mul_PValue<- round(cox4$coefficients[,5],3) 
mul_CI1<-round(cox4$conf.int[,3],3)
mul_CI2<-round(cox4$conf.int[,4],3)

mul_CI<-paste(mul_CI1,'-',mul_CI2)
mul_cox1<- data.frame("HR"=mul_HR,"CI"=mul_CI, "P"=mul_PValue)


