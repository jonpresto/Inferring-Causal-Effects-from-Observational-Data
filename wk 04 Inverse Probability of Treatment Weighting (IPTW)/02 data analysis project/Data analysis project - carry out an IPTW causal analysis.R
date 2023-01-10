library(tableone)
library(Matching)
library(MatchIt)
library(survey)

#import data
data(lalonde, package = 'MatchIt')

#create variables
age<-lalonde$age
educ<-lalonde$educ
black<-as.numeric(lalonde$race=='black')
hispan<-as.numeric(lalonde$race=='hispan')
married<-lalonde$married
nodegr<-lalonde$nodegr
re74<-lalonde$re74
re75<-lalonde$re75
re78<-lalonde$re78
treat<-lalonde$treat

#create new dataset
mydata=cbind(age,educ,black,hispan,married,nodegr,re74,re75,treat,re78)
mydata<-data.frame(mydata)

#covariates we will use
xvars<-c("age","educ","black","hispan","married","nodegr","re74","re75")

#look at a table 1
table1<- CreateTableOne(vars=xvars,strata="treat", data=mydata, test=FALSE)
## include standardized mean difference (SMD)
print(table1,smd=TRUE)

######################################
#1. propensity score model
psmodel <- glm(treat ~ age+educ+black+hispan+married+nodegr+re74+re75,
               family  = binomial(link ="logit"))

## value of propensity score for each subject
ps <-predict(psmodel, type = "response")

#create weights
weight<-ifelse(treat==1,1/(ps),1/(1-ps))

#apply weights to data
weighteddata<-svydesign(ids = ~ 1, data =mydata, weights = ~ weight)

#weighted table 1
weightedtable <-svyCreateTableOne(vars = xvars, strata = "treat", 
                                  data = weighteddata, test = FALSE)


## Show table with SMD
print(weightedtable, smd = TRUE)


#to get a weighted mean for a single covariate directly:
mean(weight[treat==1]*nodegr[treat==1])/(mean(weight[treat==1]))
mean(weight[treat==0]*nodegr[treat==0])/(mean(weight[treat==0]))


##########################################
#1. get causal risk difference
##########################################
glm.obj<-glm(re78 ~ treat, weights=weight,family=quasi(link="identity"))

#summary(glm.obj)
betaiptw<-coef(glm.obj)
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

causalrd<-(betaiptw[2])
lcl<-(betaiptw[2]-1.96*SE[2])
ucl<-(betaiptw[2]+1.96*SE[2])
c(lcl,causalrd,ucl)

# now try with truncating weights
quantile(weight, probs = c(.01, .5, .99))
th_max = quantile(weight, probs = c(.99))
th_min = quantile(weight, probs = c(.01))

truncweight<-replace(weight,weight>=th_max,th_max)
truncweight<-replace(weight,weight=<th_min,th_min)

#get causal risk difference
glm.obj<-glm(re78 ~ treat,weights=truncweight,family=quasi(link="identity"))
#summary(glm.obj)
betaiptw<-coef(glm.obj)
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

causalrd<-(betaiptw[2])
lcl<-(betaiptw[2]-1.96*SE[2])
ucl<-(betaiptw[2]+1.96*SE[2])
c(lcl,causalrd,ucl)


##########################################
#2. get causal relative risk. Weighted GLM
##########################################
glm.obj<-glm(re78 ~ treat,weights=weight,family=quasi(link=log))
#summary(glm.obj)
betaiptw<-coef(glm.obj)
#to properly account for weighting, use asymptotic (sandwich) variance
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

#get point estimate and CI for relative risk (need to exponentiate)
causalrr<-exp(betaiptw[2])
lcl<-exp(betaiptw[2]-1.96*SE[2])
ucl<-exp(betaiptw[2]+1.96*SE[2])
c(lcl,causalrr,ucl)


#######################################################
#2. fit a marginal structural model (risk difference)
msm <- (svyglm(re78~ treat, design = svydesign(~ 1, weights = ~weight, data =mydata, trunc=0.01)))
coef(msm)
confint(msm)
