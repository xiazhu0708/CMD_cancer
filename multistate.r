######Multi-state, CKB
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

dat1<-df[,c("csid","fcmd.t","fcmd.s","cmm.t","cmm.s","MN_difftime","MN_total","death_difftime","death_total",
            "smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score",
            "unhealth_score_5g","age","is_male","region_is_urban","education_2g","marital_2g",
            "income_2g","cancer_fh_2g","stroke_fh_2g","heart_fh_2g","diabetes_fh_2g",
            "cmd_fh_2g","prevalent_hyper_2g")]
dat1<-dat1 %>% mutate_at(vars("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g",
                              "unhealth_score_5g"),.funs=as.factor)
tmat<-transMat(x=list(c(2,4,5),c(3,4,5),c(4,5),c(),c()),
               names=c("Baseline","FCMD","CMM","MN","Death"))
tmat
msebmt<-msprep(data=dat1,trans=tmat,time=c(NA,"fcmd.t","cmm.t","MN_difftime","death_difftime"),
               status=c(NA,"fcmd.s","cmm.s","MN_total","death_total"),
               keep=c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score",
                      "unhealth_score_5g","age","is_male","region_is_urban","education_2g","marital_2g",
                      "income_2g","cancer_fh_2g","stroke_fh_2g","heart_fh_2g","diabetes_fh_2g",
                      "cmd_fh_2g","prevalent_hyper_2g"))
events(msebmt)

covs<-c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score",
        "unhealth_score_5g","age","is_male","region_is_urban","education_2g","marital_2g",
        "income_2g","cancer_fh_2g","stroke_fh_2g","heart_fh_2g","diabetes_fh_2g",
        "cmd_fh_2g","prevalent_hyper_2g")
msebmt<-expand.covs(msebmt,covs,longnames=FALSE)
msebmt[msebmt$id==1,]
msebmt[,c("Tstart","Tstop","time")]<-msebmt[,c("Tstart","Tstop","time")]/365.25

names(msebmt)
grep("smoking_2g.1",colnames(msebmt))#28
grep("bodyshape_2g.8",colnames(msebmt))#67
grep("unhealth_score.1",colnames(msebmt))#68
grep("unhealth_score.8",colnames(msebmt))#75
grep("age.1",colnames(msebmt))#108
grep("cancer_fh_2g.8",colnames(msebmt))#163
grep("cmd_fh_2g.1",colnames(msebmt))#188
grep("cmd_fh_2g.8",colnames(msebmt))#195
grep("prevalent_hyper_2g.8",colnames(msebmt))#203
grep("unhealth_score_5g1.1",colnames(msebmt))#76
grep("unhealth_score_5g4.8",colnames(msebmt))#107

#single LF
selected_names<-c(colnames(msebmt)[28:67],colnames(msebmt)[108:163],colnames(msebmt)[188:203])
variable_names<-paste(selected_names,collapse="+")
formula<-as.formula(paste("Surv(Tstart,Tstop,status)~",variable_names,"+strata(trans)"))
cfull<-coxph(formula,data=msebmt,method="breslow")
mul_cox<- function(x){
  mul_Beta<-round(x$coefficients[,1],3)
  mul_HR<-round(x$coefficients[,2],2)
  mul_PValue<-x$coefficients[,5]
  mul_CI1<-round(x$conf.int[,3],2)
  mul_CI2<-round(x$conf.int[,4],2)
  mul_CI<-paste0(mul_HR," (",mul_CI1,"-",mul_CI2,")")
  mul_output<-data.frame("Beta"=mul_Beta,"HR"=mul_HR,"CI5"=mul_CI1,"CI95"=mul_CI2,"HR.95CI"=mul_CI,"p"=mul_PValue)
  mul_output<-data.frame(names=row.names(mul_output),mul_output)
  colnames(mul_output)[1]<-'Characteristics'
  mul_output
}
rfmodel<-summary(cfull)
mul_cox1<-mul_cox(rfmodel)
mul_cox1

#LF score
selected_names<-c(colnames(msebmt)[76:163],colnames(msebmt)[188:203])
variable_names<-paste(selected_names,collapse="+")
formula<-as.formula(paste("Surv(Tstart,Tstop,status)~",variable_names,"+strata(trans)"))
cfull_lf<-coxph(formula,data=msebmt,method="breslow")
rfmodel_lf<-summary(cfull_lf)
mul_cox1_lf<-mul_cox(rfmodel_lf)
mul_cox1_lf

#per LF score
selected_names<-c(colnames(msebmt)[68:75],colnames(msebmt)[108:163],colnames(msebmt)[188:203])
variable_names<-paste(selected_names,collapse="+")
formula<-as.formula(paste("Surv(Tstart,Tstop,status)~",variable_names,"+strata(trans)"))
cfull_per<-coxph(formula,data=msebmt,method="breslow")
rfmodel_per<-summary(cfull_per)
mul_per<-mul_cox(rfmodel_per)
mul_per

#trans number
trans1<-subset(msebmt,trans==1)
ntrans1<-data.frame("sm"=table(trans1$status,trans1$smoking_2g)[4],"al"=table(trans1$status,trans1$alcohol_2g)[4],
                    "mt"=table(trans1$status,trans1$met_2g)[4],"dt"=table(trans1$status,trans1$diet_health_2g)[4],
                    "bs"=table(trans1$status,trans1$bodyshape_2g)[4])
lf<-table(trans1$status,trans1$unhealth_score_5g)
ntrans1<-rbind(ntrans1,lf[2,])
ntrans1$trans<-1
rownames(ntrans1)[2]<-"lf"
trans4<-subset(msebmt,trans==4)
ntrans4<-data.frame("sm"=table(trans4$status,trans4$smoking_2g)[4],"al"=table(trans4$status,trans4$alcohol_2g)[4],
                    "mt"=table(trans4$status,trans4$met_2g)[4],"dt"=table(trans4$status,trans4$diet_health_2g)[4],
                    "bs"=table(trans4$status,trans4$bodyshape_2g)[4])
lf<-table(trans4$status,trans4$unhealth_score_5g)
ntrans4<-rbind(ntrans4,lf[2,])
ntrans4$trans<-4
rownames(ntrans4)[2]<-"lf"
trans2<-subset(msebmt,trans==2)
ntrans2<-data.frame("sm"=table(trans2$status,trans2$smoking_2g)[4],"al"=table(trans2$status,trans2$alcohol_2g)[4],
                    "mt"=table(trans2$status,trans2$met_2g)[4],"dt"=table(trans2$status,trans2$diet_health_2g)[4],
                    "bs"=table(trans2$status,trans2$bodyshape_2g)[4])
lf<-table(trans2$status,trans2$unhealth_score_5g)
ntrans2<-rbind(ntrans2,lf[2,])
ntrans2$trans<-2
rownames(ntrans2)[2]<-"lf"
trans5<-subset(msebmt,trans==5)
ntrans5<-data.frame("sm"=table(trans5$status,trans5$smoking_2g)[4],"al"=table(trans5$status,trans5$alcohol_2g)[4],
                    "mt"=table(trans5$status,trans5$met_2g)[4],"dt"=table(trans5$status,trans5$diet_health_2g)[4],
                    "bs"=table(trans5$status,trans5$bodyshape_2g)[4])
lf<-table(trans5$status,trans5$unhealth_score_5g)
ntrans5<-rbind(ntrans5,lf[2,])
ntrans5$trans<-5
rownames(ntrans5)[2]<-"lf"
trans7<-subset(msebmt,trans==7)
ntrans7<-data.frame("sm"=table(trans7$status,trans7$smoking_2g)[4],"al"=table(trans7$status,trans7$alcohol_2g)[4],
                    "mt"=table(trans7$status,trans7$met_2g)[4],"dt"=table(trans7$status,trans7$diet_health_2g)[4],
                    "bs"=table(trans7$status,trans7$bodyshape_2g)[4])
lf<-table(trans7$status,trans7$unhealth_score_5g)
ntrans7<-rbind(ntrans7,lf[2,])
ntrans7$trans<-7
rownames(ntrans7)[2]<-"lf"
trans3<-subset(msebmt,trans==3)
ntrans3<-data.frame("sm"=table(trans3$status,trans3$smoking_2g)[4],"al"=table(trans3$status,trans3$alcohol_2g)[4],
                    "mt"=table(trans3$status,trans3$met_2g)[4],"dt"=table(trans3$status,trans3$diet_health_2g)[4],
                    "bs"=table(trans3$status,trans3$bodyshape_2g)[4])
lf<-table(trans3$status,trans3$unhealth_score_5g)
ntrans3<-rbind(ntrans3,lf[2,])
ntrans3$trans<-3
rownames(ntrans3)[2]<-"lf"
trans6<-subset(msebmt,trans==6)
ntrans6<-data.frame("sm"=table(trans6$status,trans6$smoking_2g)[4],"al"=table(trans6$status,trans6$alcohol_2g)[4],
                    "mt"=table(trans6$status,trans6$met_2g)[4],"dt"=table(trans6$status,trans6$diet_health_2g)[4],
                    "bs"=table(trans6$status,trans6$bodyshape_2g)[4])
lf<-table(trans6$status,trans6$unhealth_score_5g)
ntrans6<-rbind(ntrans6,lf[2,])
ntrans6$trans<-6
rownames(ntrans6)[2]<-"lf"
trans8<-subset(msebmt,trans==8)
ntrans8<-data.frame("sm"=table(trans8$status,trans8$smoking_2g)[4],"al"=table(trans8$status,trans8$alcohol_2g)[4],
                    "mt"=table(trans8$status,trans8$met_2g)[4],"dt"=table(trans8$status,trans8$diet_health_2g)[4],
                    "bs"=table(trans8$status,trans8$bodyshape_2g)[4])
lf<-table(trans8$status,trans8$unhealth_score_5g)
ntrans8<-rbind(ntrans8,lf[2,])
ntrans8$trans<-8
rownames(ntrans8)[2]<-"lf"
n.list=list(ntrans1,ntrans4,ntrans2,ntrans5,ntrans7,ntrans3,ntrans6,ntrans8)
n_com=Reduce(function(x,y){
  rbind(x,y)
},n.list)
n_com



######Multi-state, UKB
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
tmat<-transMat(x=list(c(2,4,5),c(3,4,5),c(4,5),c(),c()),
               names=c("Baseline","FCMD","CMM","MN","Death"))
tmat
msebmt<-msprep(data=dat1,trans=tmat,time=c(NA,"fcmd.t","cmm.t","MN_difftime","death_difftime"),
               status=c(NA,"fcmd.s","cmm.s","MN_total","death_total"),
               keep=c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score",
            "unhealth_score_5g","age","is_male","center","education_2g","ethnic_2g",
            "income_2g","cancer_fh_2g","stroke_fh_2g","heart_fh_2g","diabetes_fh_2g",
            "cmd_fh_2g","prevalent_hyper_2g"))
events(msebmt)

covs<-c("smoking_2g","alcohol_2g","met_2g","diet_health_2g","bodyshape_2g","unhealth_score",
        "unhealth_score_5g","age","is_male","center","education_2g","ethnic_2g",
        "income_2g","cancer_fh_2g","stroke_fh_2g","heart_fh_2g","diabetes_fh_2g",
        "cmd_fh_2g","prevalent_hyper_2g")
msebmt<-expand.covs(msebmt,covs,longnames=FALSE)
msebmt[msebmt$id==1,]
msebmt[,c("Tstart","Tstop","time")]<-msebmt[,c("Tstart","Tstop","time")]/365.25

names(msebmt)
grep("smoking_2g.1",colnames(msebmt))#28
grep("bodyshape_2g.8",colnames(msebmt))#67
grep("unhealth_score.1",colnames(msebmt))#68
grep("unhealth_score.8",colnames(msebmt))#75
grep("age.1",colnames(msebmt))#108
grep("cancer_fh_2g.8",colnames(msebmt))#195
grep("cmd_fh_2g.1",colnames(msebmt))#220
grep("cmd_fh_2g.8",colnames(msebmt))#227
grep("prevalent_hyper_2g.8",colnames(msebmt))#235
grep("unhealth_score_5g1.1",colnames(msebmt))#76
grep("unhealth_score_5g4.8",colnames(msebmt))#107
#single LF
selected_names<-c(colnames(msebmt)[28:67],colnames(msebmt)[108:195],colnames(msebmt)[220:235])
variable_names<-paste(selected_names,collapse="+")
formula<-as.formula(paste("Surv(Tstart,Tstop,status)~",variable_names,"+strata(trans)"))
cfull<-coxph(formula,data=msebmt,method="breslow")
mul_cox<- function(x){
  mul_Beta<-round(x$coefficients[,1],3)
  mul_HR<-round(x$coefficients[,2],2)
  mul_PValue<-x$coefficients[,5]
  mul_CI1<-round(x$conf.int[,3],2)
  mul_CI2<-round(x$conf.int[,4],2)
  mul_CI<-paste0(mul_HR," (",mul_CI1,"-",mul_CI2,")")
  mul_output<-data.frame("Beta"=mul_Beta,"HR"=mul_HR,"CI5"=mul_CI1,"CI95"=mul_CI2,"HR.95CI"=mul_CI,"p"=mul_PValue)
  mul_output<-data.frame(names=row.names(mul_output),mul_output)
  colnames(mul_output)[1]<-'Characteristics'
  mul_output
}
rfmodel<-summary(cfull)
mul_cox1<-mul_cox(rfmodel)
mul_cox1

#LF score
selected_names<-c(colnames(msebmt)[76:195],colnames(msebmt)[220:235])
variable_names<-paste(selected_names,collapse="+")
formula<-as.formula(paste("Surv(Tstart,Tstop,status)~",variable_names,"+strata(trans)"))
cfull_lf<-coxph(formula,data=msebmt,method="breslow")
rfmodel_lf<-summary(cfull_lf)
mul_cox1_lf<-mul_cox(rfmodel_lf)
mul_cox1_lf

#per LF score
selected_names<-c(colnames(msebmt)[68:75],colnames(msebmt)[108:195],colnames(msebmt)[220:235])
variable_names<-paste(selected_names,collapse="+")
formula<-as.formula(paste("Surv(Tstart,Tstop,status)~",variable_names,"+strata(trans)"))
cfull_per<-coxph(formula,data=msebmt,method="breslow")
rfmodel_per<-summary(cfull_per)
mul_per<-mul_cox(rfmodel_per)
mul_per

#trans number
trans1<-subset(msebmt,trans==1)
ntrans1<-data.frame("sm"=table(trans1$status,trans1$smoking_2g)[4],"al"=table(trans1$status,trans1$alcohol_2g)[4],
                    "mt"=table(trans1$status,trans1$met_2g)[4],"dt"=table(trans1$status,trans1$diet_health_2g)[4],
                    "bs"=table(trans1$status,trans1$bodyshape_2g)[4])
lf<-table(trans1$status,trans1$unhealth_score_5g)
ntrans1<-rbind(ntrans1,lf[2,])
ntrans1$trans<-1
rownames(ntrans1)[2]<-"lf"
trans4<-subset(msebmt,trans==4)
ntrans4<-data.frame("sm"=table(trans4$status,trans4$smoking_2g)[4],"al"=table(trans4$status,trans4$alcohol_2g)[4],
                    "mt"=table(trans4$status,trans4$met_2g)[4],"dt"=table(trans4$status,trans4$diet_health_2g)[4],
                    "bs"=table(trans4$status,trans4$bodyshape_2g)[4])
lf<-table(trans4$status,trans4$unhealth_score_5g)
ntrans4<-rbind(ntrans4,lf[2,])
ntrans4$trans<-4
rownames(ntrans4)[2]<-"lf"
trans2<-subset(msebmt,trans==2)
ntrans2<-data.frame("sm"=table(trans2$status,trans2$smoking_2g)[4],"al"=table(trans2$status,trans2$alcohol_2g)[4],
                    "mt"=table(trans2$status,trans2$met_2g)[4],"dt"=table(trans2$status,trans2$diet_health_2g)[4],
                    "bs"=table(trans2$status,trans2$bodyshape_2g)[4])
lf<-table(trans2$status,trans2$unhealth_score_5g)
ntrans2<-rbind(ntrans2,lf[2,])
ntrans2$trans<-2
rownames(ntrans2)[2]<-"lf"
trans5<-subset(msebmt,trans==5)
ntrans5<-data.frame("sm"=table(trans5$status,trans5$smoking_2g)[4],"al"=table(trans5$status,trans5$alcohol_2g)[4],
                    "mt"=table(trans5$status,trans5$met_2g)[4],"dt"=table(trans5$status,trans5$diet_health_2g)[4],
                    "bs"=table(trans5$status,trans5$bodyshape_2g)[4])
lf<-table(trans5$status,trans5$unhealth_score_5g)
ntrans5<-rbind(ntrans5,lf[2,])
ntrans5$trans<-5
rownames(ntrans5)[2]<-"lf"
trans7<-subset(msebmt,trans==7)
ntrans7<-data.frame("sm"=table(trans7$status,trans7$smoking_2g)[4],"al"=table(trans7$status,trans7$alcohol_2g)[4],
                    "mt"=table(trans7$status,trans7$met_2g)[4],"dt"=table(trans7$status,trans7$diet_health_2g)[4],
                    "bs"=table(trans7$status,trans7$bodyshape_2g)[4])
lf<-table(trans7$status,trans7$unhealth_score_5g)
ntrans7<-rbind(ntrans7,lf[2,])
ntrans7$trans<-7
rownames(ntrans7)[2]<-"lf"
trans3<-subset(msebmt,trans==3)
ntrans3<-data.frame("sm"=table(trans3$status,trans3$smoking_2g)[4],"al"=table(trans3$status,trans3$alcohol_2g)[4],
                    "mt"=table(trans3$status,trans3$met_2g)[4],"dt"=table(trans3$status,trans3$diet_health_2g)[4],
                    "bs"=table(trans3$status,trans3$bodyshape_2g)[4])
lf<-table(trans3$status,trans3$unhealth_score_5g)
ntrans3<-rbind(ntrans3,lf[2,])
ntrans3$trans<-3
rownames(ntrans3)[2]<-"lf"
trans6<-subset(msebmt,trans==6)
ntrans6<-data.frame("sm"=table(trans6$status,trans6$smoking_2g)[4],"al"=table(trans6$status,trans6$alcohol_2g)[4],
                    "mt"=table(trans6$status,trans6$met_2g)[4],"dt"=table(trans6$status,trans6$diet_health_2g)[4],
                    "bs"=table(trans6$status,trans6$bodyshape_2g)[4])
lf<-table(trans6$status,trans6$unhealth_score_5g)
ntrans6<-rbind(ntrans6,lf[2,])
ntrans6$trans<-6
rownames(ntrans6)[2]<-"lf"
trans8<-subset(msebmt,trans==8)
ntrans8<-data.frame("sm"=table(trans8$status,trans8$smoking_2g)[4],"al"=table(trans8$status,trans8$alcohol_2g)[4],
                    "mt"=table(trans8$status,trans8$met_2g)[4],"dt"=table(trans8$status,trans8$diet_health_2g)[4],
                    "bs"=table(trans8$status,trans8$bodyshape_2g)[4])
lf<-table(trans8$status,trans8$unhealth_score_5g)
ntrans8<-rbind(ntrans8,lf[2,])
ntrans8$trans<-8
rownames(ntrans8)[2]<-"lf"
n.list=list(ntrans1,ntrans4,ntrans2,ntrans5,ntrans7,ntrans3,ntrans6,ntrans8)
n_com=Reduce(function(x,y){
  rbind(x,y)
},n.list)
n_com



