######Cox, CKB
rm(list=ls())
options(stringsAsFactors=FALSE)
library(lubridate)
library(dplyr)
library(plyr)
library(rms)
library(survival)
library(survminer)
library(openxlsx)
library(data.table)
#load("df.rdata")

dat1<-df[,c("csid","fcmd.t","fcmd.s","cmm.t","cmm.s","MN_difftime","MN_total","death_difftime","death_total",
            "smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score",
            "unhealth_score_5g","age","is_male","region_is_urban","education_2g","marital_2g",
            "income_2g","cancer_fh_2g","stroke_fh_2g","heart_fh_2g","diabetes_fh_2g",
            "cmd_fh_2g","prevalent_hyper_2g")]
dat1<-dat1 %>% mutate_at(vars("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g",
                              "unhealth_score_5g"),.funs=as.factor)
mul_cox<- function(x){
  mul_Beta<-round(x$coefficients[,1],3)
  mul_HR<-round(x$coefficients[,2],2)
  mul_PValue<-x$coefficients[,5]
  mul_CI1<-round(x$conf.int[,3],2)
  mul_CI2<-round(x$conf.int[,4],2)
  mul_CI<-paste0(mul_HR," (",mul_CI1,"-",mul_CI2,")")
  mul_output<-data.frame("Beta"=mul_Beta,"HR"=mul_HR,"CI5"=mul_CI1,"CI95"=mul_CI2,"HR.95CI"=mul_CI,"p"=mul_PValue)
  mul_output
}

#fcmd
dat1$fcmd.etime<-dat1$fcmd.t/365.25
cox1<-coxph(Surv(fcmd.etime,fcmd.s)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
              age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
              cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox1<-summary(cox1)
mul_cox1<-mul_cox(cox1)
mul_cox1$Type<-"fcmd"
cox1_1<-coxph(Surv(fcmd.etime,fcmd.s)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g,data=dat1,x=TRUE)
cox1_1<-summary(cox1_1)
mul_cox1_1<-mul_cox(cox1_1)
mul_cox1_1$Type<-"fcmd_noHTN"
cox1lf<-coxph(Surv(fcmd.etime,fcmd.s)~unhealth_score_5g+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox1lf<-summary(cox1lf)
mul_cox1lf<-mul_cox(cox1lf)
mul_cox1lf$Type<-"lf_fcmd"
cox1lf_1<-coxph(Surv(fcmd.etime,fcmd.s)~unhealth_score_5g+
                  age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox1lf_1<-summary(cox1lf_1)
mul_cox1lf_1<-mul_cox(cox1lf_1)
mul_cox1lf_1$Type<-"lf_fcmd_noHTN"
cox1per<-coxph(Surv(fcmd.etime,fcmd.s)~unhealth_score+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox1per<-summary(cox1per)
mul_cox1per<-mul_cox(cox1per)
mul_cox1per$Type<-"per_fcmd"
cox1per_1<-coxph(Surv(fcmd.etime,fcmd.s)~unhealth_score+
                  age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox1per_1<-summary(cox1per_1)
mul_cox1per_1<-mul_cox(cox1per_1)
mul_cox1per_1$Type<-"per_fcmd_noHTN"

py<-function(x){
  event<-x$event
  pyears<-x$pyears
  df<-data.frame(event,pyears)
  df$py<-df$event/df$pyears*100000
  df
}
variables<-c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score_5g")
result_list<-list()
for (var in variables) {
  formula<-as.formula(paste("Surv(fcmd.etime,fcmd.s)~",var))
  fit<-pyears(formula,data=dat1,scale=1)
  pytable<-py(fit)
  pytable$variable<-var
  result_list[[var]]<-pytable
}
fcmd_result<-do.call(rbind,result_list)
fcmd_result$type<-"fcmd"

#cmm
dat1$cmm.etime<-dat1$cmm.t/365.25
cox2<-coxph(Surv(cmm.etime,cmm.s)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
              age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
              cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox2<-summary(cox2)
mul_cox2<-mul_cox(cox2)
mul_cox2$Type<-"cmm"
cox2_1<-coxph(Surv(cmm.etime,cmm.s)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g,data=dat1,x=TRUE)
cox2_1<-summary(cox2_1)
mul_cox2_1<-mul_cox(cox2_1)
mul_cox2_1$Type<-"cmm_noHTN"
cox2lf<-coxph(Surv(cmm.etime,cmm.s)~unhealth_score_5g+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox2lf<-summary(cox2lf)
mul_cox2lf<-mul_cox(cox2lf)
mul_cox2lf$Type<-"lf_cmm"
cox2lf_1<-coxph(Surv(cmm.etime,cmm.s)~unhealth_score_5g+
                  age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox2lf_1<-summary(cox2lf_1)
mul_cox2lf_1<-mul_cox(cox2lf_1)
mul_cox2lf_1$Type<-"lf_cmm_noHTN"
cox2per<-coxph(Surv(cmm.etime,cmm.s)~unhealth_score+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox2per<-summary(cox2per)
mul_cox2per<-mul_cox(cox2per)
mul_cox2per$Type<-"per_cmm"
cox2per_1<-coxph(Surv(cmm.etime,cmm.s)~unhealth_score+
                  age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox2per_1<-summary(cox2per_1)
mul_cox2per_1<-mul_cox(cox2per_1)
mul_cox2per_1$Type<-"per_cmm_noHTN"

variables<-c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score_5g")
result_list<-list()
for (var in variables) {
  formula<-as.formula(paste("Surv(cmm.etime,cmm.s)~",var))
  fit<-pyears(formula,data=dat1,scale=1)
  pytable<-py(fit)
  pytable$variable<-var
  result_list[[var]]<-pytable
}
cmm_result<-do.call(rbind,result_list)
cmm_result$type<-"cmm"

#mn
dat1$MN_etime<-dat1$MN_difftime/365.25
cox3<-coxph(Surv(MN_etime,MN_total)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
              age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
              cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox3<-summary(cox3)
mul_cox3<-mul_cox(cox3)
mul_cox3$Type<-"mn"
cox3_1<-coxph(Surv(MN_etime,MN_total)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g,data=dat1,x=TRUE)
cox3_1<-summary(cox3_1)
mul_cox3_1<-mul_cox(cox3_1)
mul_cox3_1$Type<-"mn_noHTN"
cox3lf<-coxph(Surv(MN_etime,MN_total)~unhealth_score_5g+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox3lf<-summary(cox3lf)
mul_cox3lf<-mul_cox(cox3lf)
mul_cox3lf$Type<-"lf_mn"
cox3lf_1<-coxph(Surv(MN_etime,MN_total)~unhealth_score_5g+
                  age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox3lf_1<-summary(cox3lf_1)
mul_cox3lf_1<-mul_cox(cox3lf_1)
mul_cox3lf_1$Type<-"lf_mn_noHTN"
cox3per<-coxph(Surv(MN_etime,MN_total)~unhealth_score+
                age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox3per<-summary(cox3per)
mul_cox3per<-mul_cox(cox3per)
mul_cox3per$Type<-"per_mn"
cox3per_1<-coxph(Surv(MN_etime,MN_total)~unhealth_score+
                  age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox3per_1<-summary(cox3per_1)
mul_cox3per_1<-mul_cox(cox3per_1)
mul_cox3per_1$Type<-"per_mn_noHTN"

variables<-c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score_5g")
result_list<-list()
for (var in variables) {
  formula<-as.formula(paste("Surv(MN_etime,MN_total)~",var))
  fit<-pyears(formula,data=dat1,scale=1)
  pytable<-py(fit)
  pytable$variable<-var
  result_list[[var]]<-pytable
}
mn_result<-do.call(rbind,result_list)
mn_result$type<-"mn"

cox.list=list(mul_cox1,mul_cox1_1,mul_cox1lf,mul_cox1lf_1,mul_cox1per,mul_cox1per_1,
              mul_cox2,mul_cox2_1,mul_cox2lf,mul_cox2lf_1,mul_cox2per,mul_cox2per_1,
              mul_cox3,mul_cox3_1,mul_cox3lf,mul_cox3lf_1,mul_cox3per,mul_cox3per_1)
cox_com=Reduce(function(x,y){
  rbind(x,y)
},cox.list)
cox_com

py.list=list(fcmd_result,cmm_result,mn_result)
py_com=Reduce(function(x,y){
  rbind(x,y)
},py.list)
py_com

#time dependent fcmd&cmm
dat1$fcmd.t2<-ifelse(dat1$fcmd.s==0,NA,dat1$fcmd.t)
dat1$cmm.t2<-ifelse(dat1$cmm.s==0,NA,dat1$cmm.t)
summary(dat1$fcmd.t2)
table(dat1$fcmd.s)
summary(dat1$cmm.t2)
table(dat1$cmm.s)
names(dat1)
newdata<-tmerge(data1=dat1[,c(1,10:28)],data2=dat1,id=csid,MN=event(MN_difftime,MN_total))
newdata<-tmerge(newdata,dat1,id=csid,fcmd=tdc(fcmd.t2))
newdata<-tmerge(newdata,dat1,id=csid,cmm=tdc(cmm.t2))
newdata<-tmerge(newdata,newdata,id=csid,cmdnum=cumtdc(tstart))
attr(newdata,"tcount")
names(newdata)
table(newdata$fcmd,exclude=NULL)
table(newdata$cmm,exclude=NULL)
table(newdata$cmdnum,exclude=NULL)
newdata$cmd<-ifelse(newdata$fcmd==0&newdata$cmm==0,0,
             ifelse(newdata$fcmd==1&newdata$cmm==0,1,
             ifelse(newdata$fcmd==1&newdata$cmm==1,2,NA)))
table(newdata$cmd,exclude=NULL)
newdata$cmd<-as.factor(newdata$cmd)
py<-function(x){
  event<-x$event
  pyears<-x$pyears
  df<-data.frame(event,pyears)
  df$py<-df$event/df$pyears*100000
  df
}
fit0<-pyears(Surv(tstart/365.25,tstop/365.25,MN)~cmd,data=newdata,scale=1)
py(fit0)

tcox1<-coxph(Surv(tstart,tstop,MN)~cmd,data=newdata,ties="breslow")
tcox1<-summary(tcox1)
mul_tcox1<-mul_cox(tcox1)
mul_tcox1$Type<-"cmd_unadj"
tcox1_1<-coxph(Surv(tstart,tstop,MN)~cmd+smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
                 age+is_male+region_is_urban+education_2g+marital_2g+income_2g+cancer_fh_2g+
                 cmd_fh_2g+prevalent_hyper_2g,data=newdata,ties="breslow")
tcox1_1<-summary(tcox1_1)
mul_tcox1_1<-mul_cox(tcox1_1)
mul_tcox1_1$Type<-"cmd_adj"

cox.list=list(mul_tcox1,mul_tcox1_1)
cox_com=Reduce(function(x,y){
  rbind(x,y)
},cox.list)
cox_com



######Cox, UKB
rm(list=ls())
options(stringsAsFactors=FALSE)
library(mstate)
library(lubridate)
library(dplyr)
library(plyr)
library(rms)
library(survival)
library(survminer)
library(openxlsx)
library(data.table)
#load("df.rdata")

dat1<-df[,c("f.eid","fcmd.t","fcmd.s","cmm.t","cmm.s","MN_difftime","MN_total","death_difftime","death_total",
            "smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score",
            "unhealth_score_5g","age","is_male","center","education_2g","ethnic_2g",
            "income_2g","cancer_fh_2g","stroke_fh_2g","heart_fh_2g","diabetes_fh_2g",
            "cmd_fh_2g","prevalent_hyper_2g")]
dat1<-dat1 %>% mutate_at(vars("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g",
            "unhealth_score_5g","is_male","center","education_2g","ethnic_2g",
            "income_2g","cancer_fh_2g","stroke_fh_2g","heart_fh_2g","diabetes_fh_2g",
            "cmd_fh_2g","prevalent_hyper_2g"),.funs=as.factor)
mul_cox<- function(x){
  mul_Beta<-round(x$coefficients[,1],3)
  mul_HR<-round(x$coefficients[,2],2)
  mul_PValue<-x$coefficients[,5]
  mul_CI1<-round(x$conf.int[,3],2)
  mul_CI2<-round(x$conf.int[,4],2)
  mul_CI<-paste0(mul_HR," (",mul_CI1,"-",mul_CI2,")")
  mul_output<-data.frame("Beta"=mul_Beta,"HR"=mul_HR,"CI5"=mul_CI1,"CI95"=mul_CI2,"HR.95CI"=mul_CI,"p"=mul_PValue)
  mul_output
}

#fcmd
dat1$fcmd.etime<-dat1$fcmd.t/365.25
cox1<-coxph(Surv(fcmd.etime,fcmd.s)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
              age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
              cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox1<-summary(cox1)
mul_cox1<-mul_cox(cox1)
mul_cox1$Type<-"fcmd"
cox1_1<-coxph(Surv(fcmd.etime,fcmd.s)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g,data=dat1,x=TRUE)
cox1_1<-summary(cox1_1)
mul_cox1_1<-mul_cox(cox1_1)
mul_cox1_1$Type<-"fcmd_noHTN"
cox1lf<-coxph(Surv(fcmd.etime,fcmd.s)~unhealth_score_5g+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox1lf<-summary(cox1lf)
mul_cox1lf<-mul_cox(cox1lf)
mul_cox1lf$Type<-"lf_fcmd"
cox1lf_1<-coxph(Surv(fcmd.etime,fcmd.s)~unhealth_score_5g+
                  age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox1lf_1<-summary(cox1lf_1)
mul_cox1lf_1<-mul_cox(cox1lf_1)
mul_cox1lf_1$Type<-"lf_fcmd_noHTN"
cox1per<-coxph(Surv(fcmd.etime,fcmd.s)~unhealth_score+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox1per<-summary(cox1per)
mul_cox1per<-mul_cox(cox1per)
mul_cox1per$Type<-"per_fcmd"
cox1per_1<-coxph(Surv(fcmd.etime,fcmd.s)~unhealth_score+
                  age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox1per_1<-summary(cox1per_1)
mul_cox1per_1<-mul_cox(cox1per_1)
mul_cox1per_1$Type<-"per_fcmd_noHTN"

py<-function(x){
  event<-x$event
  pyears<-x$pyears
  df<-data.frame(event,pyears)
  df$py<-df$event/df$pyears*100000
  df
}
variables<-c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score_5g")
result_list<-list()
for (var in variables) {
  formula<-as.formula(paste("Surv(fcmd.etime,fcmd.s)~",var))
  fit<-pyears(formula,data=dat1,scale=1)
  pytable<-py(fit)
  pytable$variable<-var
  result_list[[var]]<-pytable
}
fcmd_result<-do.call(rbind,result_list)
fcmd_result$type<-"fcmd"

#cmm
dat1$cmm.etime<-dat1$cmm.t/365.25
cox2<-coxph(Surv(cmm.etime,cmm.s)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
              age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
              cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox2<-summary(cox2)
mul_cox2<-mul_cox(cox2)
mul_cox2$Type<-"cmm"
cox2_1<-coxph(Surv(cmm.etime,cmm.s)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g,data=dat1,x=TRUE)
cox2_1<-summary(cox2_1)
mul_cox2_1<-mul_cox(cox2_1)
mul_cox2_1$Type<-"cmm_noHTN"
cox2lf<-coxph(Surv(cmm.etime,cmm.s)~unhealth_score_5g+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox2lf<-summary(cox2lf)
mul_cox2lf<-mul_cox(cox2lf)
mul_cox2lf$Type<-"lf_cmm"
cox2lf_1<-coxph(Surv(cmm.etime,cmm.s)~unhealth_score_5g+
                  age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox2lf_1<-summary(cox2lf_1)
mul_cox2lf_1<-mul_cox(cox2lf_1)
mul_cox2lf_1$Type<-"lf_cmm_noHTN"
cox2per<-coxph(Surv(cmm.etime,cmm.s)~unhealth_score+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox2per<-summary(cox2per)
mul_cox2per<-mul_cox(cox2per)
mul_cox2per$Type<-"per_cmm"
cox2per_1<-coxph(Surv(cmm.etime,cmm.s)~unhealth_score+
                  age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox2per_1<-summary(cox2per_1)
mul_cox2per_1<-mul_cox(cox2per_1)
mul_cox2per_1$Type<-"per_cmm_noHTN"

variables<-c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score_5g")
result_list<-list()
for (var in variables) {
  formula<-as.formula(paste("Surv(cmm.etime,cmm.s)~",var))
  fit<-pyears(formula,data=dat1,scale=1)
  pytable<-py(fit)
  pytable$variable<-var
  result_list[[var]]<-pytable
}
cmm_result<-do.call(rbind,result_list)
cmm_result$type<-"cmm"

#mn
dat1$MN_etime<-dat1$MN_difftime/365.25
cox3<-coxph(Surv(MN_etime,MN_total)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
              age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
              cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox3<-summary(cox3)
mul_cox3<-mul_cox(cox3)
mul_cox3$Type<-"mn"
cox3_1<-coxph(Surv(MN_etime,MN_total)~smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g,data=dat1,x=TRUE)
cox3_1<-summary(cox3_1)
mul_cox3_1<-mul_cox(cox3_1)
mul_cox3_1$Type<-"mn_noHTN"
cox3lf<-coxph(Surv(MN_etime,MN_total)~unhealth_score_5g+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox3lf<-summary(cox3lf)
mul_cox3lf<-mul_cox(cox3lf)
mul_cox3lf$Type<-"lf_mn"
cox3lf_1<-coxph(Surv(MN_etime,MN_total)~unhealth_score_5g+
                  age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox3lf_1<-summary(cox3lf_1)
mul_cox3lf_1<-mul_cox(cox3lf_1)
mul_cox3lf_1$Type<-"lf_mn_noHTN"
cox3per<-coxph(Surv(MN_etime,MN_total)~unhealth_score+
                age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                cmd_fh_2g+prevalent_hyper_2g,data=dat1,x=TRUE)
cox3per<-summary(cox3per)
mul_cox3per<-mul_cox(cox3per)
mul_cox3per$Type<-"per_mn"
cox3per_1<-coxph(Surv(MN_etime,MN_total)~unhealth_score+
                  age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                  cmd_fh_2g,data=dat1,x=TRUE)
cox3per_1<-summary(cox3per_1)
mul_cox3per_1<-mul_cox(cox3per_1)
mul_cox3per_1$Type<-"per_mn_noHTN"

variables<-c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score_5g")
result_list<-list()
for (var in variables) {
  formula<-as.formula(paste("Surv(MN_etime,MN_total)~",var))
  fit<-pyears(formula,data=dat1,scale=1)
  pytable<-py(fit)
  pytable$variable<-var
  result_list[[var]]<-pytable
}
mn_result<-do.call(rbind,result_list)
mn_result$type<-"mn"

cox.list=list(mul_cox1,mul_cox1_1,mul_cox1lf,mul_cox1lf_1,mul_cox1per,mul_cox1per_1,
              mul_cox2,mul_cox2_1,mul_cox2lf,mul_cox2lf_1,mul_cox2per,mul_cox2per_1,
              mul_cox3,mul_cox3_1,mul_cox3lf,mul_cox3lf_1,mul_cox3per,mul_cox3per_1)
cox_com=Reduce(function(x,y){
  rbind(x,y)
},cox.list)
cox_com

py.list=list(fcmd_result,cmm_result,mn_result)
py_com=Reduce(function(x,y){
  rbind(x,y)
},py.list)
py_com

#time dependent fcmd&cmm
dat1$fcmd.t2<-ifelse(dat1$fcmd.s==0,NA,dat1$fcmd.t)
dat1$cmm.t2<-ifelse(dat1$cmm.s==0,NA,dat1$cmm.t)
summary(dat1$fcmd.t2)
table(dat1$fcmd.s)
summary(dat1$cmm.t2)
table(dat1$cmm.s)
names(dat1)
newdata<-tmerge(data1=dat1[,c(1,10:28)],data2=dat1,id=f.eid,MN=event(MN_difftime,MN_total))
newdata<-tmerge(newdata,dat1,id=f.eid,fcmd=tdc(fcmd.t2))
newdata<-tmerge(newdata,dat1,id=f.eid,cmm=tdc(cmm.t2))
newdata<-tmerge(newdata,newdata,id=f.eid,cmdnum=cumtdc(tstart))
attr(newdata,"tcount")
names(newdata)
table(newdata$fcmd,exclude=NULL)
table(newdata$cmm,exclude=NULL)
table(newdata$cmdnum,exclude=NULL)
newdata$cmd<-ifelse(newdata$fcmd==0&newdata$cmm==0,0,
             ifelse(newdata$fcmd==1&newdata$cmm==0,1,
             ifelse(newdata$fcmd==1&newdata$cmm==1,2,NA)))
table(newdata$cmd,exclude=NULL)
newdata$cmd<-as.factor(newdata$cmd)
py<-function(x){
  event<-x$event
  pyears<-x$pyears
  df<-data.frame(event,pyears)
  df$py<-df$event/df$pyears*100000
  df
}
fit0<-pyears(Surv(tstart/365.25,tstop/365.25,MN)~cmd,data=newdata,scale=1)
py(fit0)

tcox1<-coxph(Surv(tstart,tstop,MN)~cmd,data=newdata,ties="breslow")
tcox1<-summary(tcox1)
mul_tcox1<-mul_cox(tcox1)
mul_tcox1$Type<-"cmd_unadj"
tcox1_1<-coxph(Surv(tstart,tstop,MN)~cmd+smoking_2g+alcohol_2g+met_2g+diet_health_2g+bodyshape_2g+
                 age+is_male+center+education_2g+ethnic_2g+income_2g+cancer_fh_2g+
                 cmd_fh_2g+prevalent_hyper_2g,data=newdata,ties="breslow")
tcox1_1<-summary(tcox1_1)
mul_tcox1_1<-mul_cox(tcox1_1)
mul_tcox1_1$Type<-"cmd_adj"

cox.list=list(mul_tcox1,mul_tcox1_1)
cox_com=Reduce(function(x,y){
  rbind(x,y)
},cox.list)
cox_com
